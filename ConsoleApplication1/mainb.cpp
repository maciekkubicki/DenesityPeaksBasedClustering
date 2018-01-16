// ConsoleApplication1.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <flann/mpi/queries.h>
#include <flann/mpi/index.h>
#include <boost/thread/thread.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/regex.hpp>

#include <flann/flann.hpp>
#include <fstream>
#include <sstream>
#include <limits>
#include <iterator>
#include <iostream>
#include <algorithm>

enum outputColumn
{
	R = 0,//Global denesity
	D = 1,//Minimum Distance to neighbour with bigger rho
	N = 2,//Nearest neighbour with bigger rho
	C = 3,//Cluster number - default 0
	CR = 4,//Local denesity(only same cluster neighbours)
	iB = 5,//is border region
};


using namespace std;
using namespace flann;

struct clusterInfo
{
	int clusterIndex;
	int maxDenesity;
};

typedef struct clusterInfo cInfo;


template <typename T>
void zeroMatrix(Matrix<T> mat, int rows, int cols)
{
	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < cols; ++j)
			mat[i][j] = T(0);
}


void setMinMax(Matrix<float> output, int rows, float *min, float *max)
{
	*min = std::numeric_limits<float>::max();
	*max = std::numeric_limits<float>::min();
	for (int i = 0; i < rows; ++i)
	{
		if (output[i][D] < *min)
			*min = output[i][D];
		if (output[i][D] > *max)
			*max = output[i][D];
	}
}

int pickClusterCenters(Matrix<float> output, int rows, int cols, float minDist, int minDens)
{
	std::cout << "Picing cluster centers" << std::endl;
	int number = 1;
	for (int i = 0; i < rows; ++i)
		if (output[i][R] > minDens && output[i][D] > minDist || output[i][N] == INT_MAX)
		{
			cout << "IDX " << i << endl;
			output[i][C] = number++;
		}
	return --number;
}

int getClusterNumberFromNeighbour(Matrix<float> output, int neighbour)
{
	//cout << neighbour << " ";
	if (output[neighbour][C] != 0)
	{
		//cout << "RET " << output[neighbour][C];
		return output[neighbour][C];
	}
	else return getClusterNumberFromNeighbour(output, output[neighbour][N]);
}
void pickClusters(Matrix<float> output, int rows, int cols)
{
	std::cout << "Picking clusters" << std::endl;
	bool hasCenters = false;
	for (int i = 0; i < rows; ++i)
		if (output[i][C] != 0)
			hasCenters = true;
	if (!hasCenters)
	{
		std::cout << "Bad parameters - no cluster centers picked" << std::endl;
		return;
	}
	for (int i = 0; i < rows; ++i)
	{
		//cout << "pkt " << i << " ";
		if (output[i][C] == 0)
			output[i][C] = getClusterNumberFromNeighbour(output, output[i][N]);
		//cout << endl;
	}
}

void normalizeDist(Matrix<float> output, int rows, int cols, float min, float max)
{
	std::cout << "Normalize distances" << std::endl;
	for (int i = 0; i < rows; ++i)
		output[i][D] = (output[i][D] - min) / (max - min);
}


float findMaxDist(Matrix<float> data, int rows, int cols)
{
	float max = INT_MIN;
	float s = 0.f;
	for (int i = 0; i < rows; ++i)
	{
		for (int j = i + 1; j < rows; ++j)
		{
			s = 0.f;
			for (int k = 0; k < cols; ++k)
				s += pow(data[j][k] - data[i][k], 2);
			if (s > max)
				max = s;
		}
	}
	return sqrt(max);//sqrt powinno byc
}

float findMinDist(Matrix<float> data, int rows, int cols)
{
	float min = INT_MAX;
	float s = 0.f;
	for (int i = 0; i < rows; ++i)
	{
		for (int j = i + 1; j < rows; ++j)
		{
			s = 0.f;
			for (int k = 0; k < cols; ++k)
				s += pow(data[j][k] - data[i][k], 2);
			if (s < min)
				min = s;
		}
	}
	return sqrt(min);//sqrt powinno byc
}

//funkcja dzieli macierz na dane i etykiety(ostatnia kolumna to etykiety)
Matrix<float> splitData(Matrix<float> data, Matrix<int> labels, int rows, int cols)
{
	Matrix<float> newData = Matrix<float>(new float[rows*(cols - 1)], rows, cols - 1);

	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols - 1; ++j)
			newData[i][j] = data[i][j];
		labels[i][0] = (int)data[i][cols - 1];
	}

	return newData;
}


template <typename T>
void writeData(const char *fname, Matrix<T> mat, int rows, int cols)
{
	std::cout << "Writing data" << std::endl;
	ofstream fp(fname, std::ofstream::out);
	if (fp.is_open())
	{
		for (int i = 0; i < rows; ++i)
		{
			stringstream strline(stringstream::in | stringstream::out);
			for (int j = 0; j < cols; ++j)
				strline << mat[i][j] << " ";
			if (i < rows - 1)
				strline << std::endl;
			fp << strline.str();
			strline.clear();
		}
		fp.close();
	}
}
float maxx(Matrix<float> dists, int rows, int maxR)
{
	float m = std::numeric_limits<float>::min();
	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < maxR; ++j)
			if (dists[i][j] > m)
				m = dists[i][j];
	return m;
}
void getMinDist(Matrix<float> mat, Matrix<float> output, int rows, int cols, int maxR)
{
	std::cout << "Getting min distance from higher denesity" << std::endl;
	float min;
	float max = findMaxDist(mat, rows, cols);
	//setMinMax(output, rows, &min, &max);
	Matrix<float> q(new float[cols], 1, cols);
	Index<L2<float> > index(mat, flann::KDTreeIndexParams(4));
	index.buildIndex();
	Matrix<int> indices = Matrix<int>(new int[maxR*rows], rows, maxR);
	Matrix<float> dists = Matrix<float>(new float[maxR*rows], rows, maxR);
	index.knnSearch(mat, indices, dists, maxR, flann::SearchParams(128));
	max = maxx(dists, rows, maxR);
	for (int i = 0; i < rows; ++i)
	{
		if (i % 500 == 0)
			cout << i << endl;
		//for (int j = 0; j < cols; ++j)
		//	q[0][j] = mat[i][j];
		//Matrix<int> indices = Matrix<int>(new int[maxR], 1, maxR);
		//Matrix<float> dists = Matrix<float>(new float[maxR], 1, maxR);
		//index.knnSearch(q, indices, dists, maxR, flann::SearchParams(128));

		for (int k = 1; k < maxR; ++k)
		{
			//cout << i << " " << k << " " << indices[0][k] << " " << dists[0][k] << endl;
			int nn = indices[i][k];
			if (output[i][R] < output[nn][R])
			{
				output[i][D] = dists[i][k];
				output[i][N] = nn;
				break;
			}
			else if (k == maxR - 1)
			{

				output[i][D] = max;
				output[i][N] = INT_MAX;
				//cout << "Alert " << endl;
				//cout << i << " " << output[i][R] << " " << k << " " << indices[i][k] << " " << dists[i][k] << endl;

			}
		}




	}
	delete[] indices.ptr();
	delete[] dists.ptr();
	delete[] q.ptr();
}

int countMaxes(Matrix<float> output, int rows, int *idx = NULL)
{

	int i = 0;
	for (int j = 0; j < rows; ++j)
		if (output[j][N] == INT_MAX)
		{
			++i;
			if (idx != NULL)
				*idx = j;
		}
	return i;
}



void fix2(Matrix<float> mat, Matrix<float> output, float radius, int rows, int cols, int maxR)
{
	std::cout << "Fixing " << std::endl;
	bool changed = false;
	Index<L2<float> > index(mat, flann::KDTreeIndexParams(4));
	index.buildIndex();
	Matrix<float> q(new float[cols], 1, cols);
	int maxRho = INT_MIN;
	int _num = 0;
	int _num2 = 0;
	int ix = 0;
	cout << "Przed " << endl;
	for (int i = 0; i < rows; ++i)
	{
		if (output[i][R] > maxRho)
		{
			ix = i;
			_num = 1;
			maxRho = output[i][R];
		}
		else if (output[i][R] == maxRho)
			++_num;

		if (output[i][N] == INT_MAX)
			++_num2;
	}
	cout << "num " << _num << " " << _num2 << endl;
	if (_num2 < 1)
		return;
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
			q[0][j] = mat[i][j];
		if (output[i][N] == INT_MAX && (i != ix && _num != 1))
		{
			cout << "Tak" << endl;
			Matrix<int> indices = Matrix<int>(new int[maxR], 1, maxR);
			Matrix<float> dists = Matrix<float>(new float[maxR], 1, maxR);
			zeroMatrix(indices, 1, maxR);
			zeroMatrix(dists, 1, maxR);
			int num = 0;
			int nn = index.radiusSearch(q, indices, dists, radius, flann::SearchParams(128));
			for (int k = 0; k < nn; ++k)
			{
				if (output[indices[0][k]][R] == output[i][R])
					++num;
				else if (output[indices[0][k]][R] > output[i][R])
				{
					num = 0;
					break;
				}
			}
			if (num > 1)
			{
				cout << "Zmiana" << endl;
				++output[i][R];
				if (output[i][R] == maxRho)
					++output[ix][R];
				changed = true;
			}
		}
	}

	if (changed)
		getMinDist(mat, output, rows, cols, maxR);
	int m = countMaxes(output, rows, &ix);
	if (m > 1)
	{
		float min, max;
		setMinMax(output, rows, &min, &max);
		++output[ix][R];
		output[ix][D] = max;
		getMinDist(mat, output, rows, cols, maxR);
	}
}



void getDenesity(Matrix<float> mat, Matrix<float> output, float radius, int rows, int cols, int maxR)
{
	std::cout << "Getting denesity" << std::endl;
	Matrix<float> q(new float[cols], 1, cols);
	Index<L2<float> > index(mat, flann::LinearIndexParams());
	index.buildIndex();
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
			q[0][j] = mat[i][j];

		Matrix<int> indices = Matrix<int>(new int[maxR], 1, maxR);
		Matrix<float> dists = Matrix<float>(new float[maxR], 1, maxR);
		output[i][R] = index.radiusSearch(q, indices, dists, radius, flann::SearchParams(128));
		delete[] indices.ptr();
		delete[] dists.ptr();
		output[i][3] = 0;


	}
	delete[] q.ptr();
}

Matrix<float> readData(const char *fname, int *rows, int *cols)
{
	std::cout << "Reading data" << std::endl;
	Matrix<float> mat;
	*rows = 0;
	*cols = 0;
	float test[100];
	ifstream fp(fname);
	if (fp.is_open())
	{
		string line;
		stringstream strline(stringstream::in | stringstream::out);

		while (!fp.eof())
		{
			getline(fp, line);
			//cout << line << endl;

			if (*rows == 0)
			{
				strline << line;
				while (strline >> test[*cols])
				{
					//cout << test[*cols] << endl;
					++*cols;
				}
			}
			strline.clear();
			++*rows;
		}
		mat = Matrix<float>(new float[*rows * *cols], *rows, *cols);
		fp.clear();
		fp.seekg(0, ios::beg);
		int i = 0, j = 0;
		while (!fp.eof())
		{
			getline(fp, line);

			strline << line;
			j = 0;
			while (strline >> mat[i][j])
			{
				//cout << mat[i][j] << " ";
				++j;
			}
			//cout << endl;
			strline.clear();
			++i;
		}

		fp.close();
	}
	std::cout << "RD " << *rows << " " << *cols << std::endl;
	return mat;
}

void distanceToOtherCluster(Matrix<float> mat, Matrix<float> output, int rows, int cols, int maxR, float radius)
{
	std::cout << "Counting distance to other cluster" << std::endl;
	Matrix<float> q(new float[cols], 1, cols);
	Index<L2<float> > index(mat, flann::KDTreeIndexParams(4));
	index.buildIndex();
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
			q[0][j] = mat[i][j];

		Matrix<int> indices = Matrix<int>(new int[maxR], 1, maxR);
		Matrix<float> dists = Matrix<float>(new float[maxR], 1, maxR);
		int nn = index.radiusSearch(q, indices, dists, radius, flann::SearchParams(128));
		int ix = output[i][C];
		output[i][CR] = 0;
		output[i][iB] = 0;
		for (int k = 0; k < nn; ++k)
		{
			int ix2 = output[indices[0][k]][C];
			if (ix != ix2 && !output[i][iB])
			{
				//std::cout << "Border point nr " << i << " bec " << ix << " and " << ix2 << std::endl;
				output[i][iB] = 1;
			}
			++output[i][CR];

		}

	}
}

std::vector<clusterInfo> getMaxBorderDenesity(Matrix<float> output, int numberOfClusters, int rows, outputColumn E)
{
	if (E != CR && E != R)
		E = R;
	std::vector<int> dens(numberOfClusters + 1, 0);
	for (int r = 0; r < rows; ++r)
	{
		for (int i = 1; i <= numberOfClusters; ++i)
		{
			if (output[r][C] == i && dens[i] < output[r][E] && output[r][iB])
			{
				dens[i] = output[r][E];
				break;
			}
		}
	}
	std::vector<clusterInfo> c;
	for (int i = 0; i < numberOfClusters; ++i)
	{
		clusterInfo ci{ i + 1,dens[i + 1] };
		c.push_back(ci);
	}
	return c;
}

void excludeHaloPoints(Matrix<float> output, std::vector<clusterInfo> cInfo, int numberOfClusters, int rows, float factor = 1.f)
{
	std::cout << "Excluding halo points" << std::endl;
	for (int r = 0; r < rows; ++r)
	{
		for (int i = 1; i <= numberOfClusters; ++i)
		{
			if (output[r][C] == i && output[r][R] < cInfo[i - 1].maxDenesity / factor)
				output[r][C] = 0;
		}
	}
}
std::vector<int> countRowsPerNode(int size, int rows)
{
	std::vector<int> rpn(size, 0);
	int i = 0;
	while (i < rows)
	{
		++rpn[i%size];
		++i;
	}
	return rpn;
}

void clustrering(const char *fileName, float radius, float minDistF, int minDensF, bool excludeHalo = false, outputColumn densCol = R, float factor = 1.f)
{
	int nn = 5;
	int a = 2;
	int maxR = 10000;
	int rows, cols;
	//float radius = 4.2;//0.035//4.5//0.9
	//dataset2 30 0.1 1%
	//iris 0.5, 0.25, 30// 
	//seeds_dataset 0.9, 0.2, 15// 0.15 10
	//myFile 
	//gener 5, 0.02, 50
	//gener2 1.5 0.015, 200 - ciekawy
	//gener2 4.2, 0.02, 250
	//gener2 2.4, 0.02, 250
	//gener2 2.3, 0.01, 250

	Matrix<float> dataset = readData("gener2.txt", &rows, &cols);
	maxR = rows;
	double start = clock();
	Matrix<float> extraData(new float[rows * 6], rows, 6);
	float max, min;
	//max = findMaxDist(dataset, rows, cols);
	//min = findMinDist(dataset, rows, cols);

	//cout << "Euc Max " << max << " Min " << min << endl;

	//float minDistF = 0.02;
	//int minDensF = 250;//rows / 100 * 10;
	getDenesity(dataset, extraData, radius, rows, cols, rows);
	getMinDist(dataset, extraData, rows, cols, maxR);
	while (countMaxes(extraData, rows) != 1)
	{
		cout << countMaxes(extraData, rows);

		fix2(dataset, extraData, radius, rows, cols, maxR);
	}
	setMinMax(extraData, rows, &min, &max);
	cout << "Got Max " << max << " Min " << min << endl;
	normalizeDist(extraData, rows, cols, min, max);
	int numOfClusters = pickClusterCenters(extraData, rows, cols, minDistF, minDensF);
	cout << "NUM " << numOfClusters << endl;
	pickClusters(extraData, rows, cols);
	if (excludeHalo)
	{
		distanceToOtherCluster(dataset, extraData, rows, cols, rows, radius);
		std::vector<clusterInfo> ci = getMaxBorderDenesity(extraData, numOfClusters, rows, densCol);
		for (int i = 0; i < numOfClusters; ++i)
			std::cout << ci[i].clusterIndex << " " << ci[i].maxDenesity << std::endl;
		excludeHaloPoints(extraData, ci, numOfClusters, rows, factor);
	}
	writeData<float>("output.dat", extraData, rows, 6);
	cout << rows << " " << cols << endl;
	cout << "NUM " << numOfClusters << endl;
	//cout << "NUM " << numOfClusters;
	//Matrix<int> labels(new int[rows*1],rows,1);
	//Matrix<float> data = splitData(dataset, labels, rows, cols);
	//--cols;
	//cout << "tu" << endl;
	//writeData<float>("data.dat", data, rows, cols);
	//writeData<int>("labels.dat", labels, rows, 1);
	//delete[] labels.ptr();
	//delete[] data.ptr();
	delete[] dataset.ptr();
	delete[] extraData.ptr();
	//delete[] q.ptr();
	double end = clock();
	double diff = (double)(end - start) / CLOCKS_PER_SEC;
	std::cout << "Time " << diff << std::endl;
}

int main(int argc, char **argv)
{
	const char *fname = NULL;
	float radius = 0.f;
	float minDist = 1.f;
	int minDens = 10;
	bool exHalo = false;
	outputColumn densType = R;
	float factor = 1.f;
	std::cout << "Num of params " << argc << std::endl;
	if (argc < 5)
	{
		std::cout << "Input at least 4 first arguments of  filename(const char*), radius(float), minDistance(0-1 float), minimumDenesity(int), excludeHalo(bool), densType(0/4), excludingFactor(>=1)" << std::endl;
		exit(-1);
	}
	fname = argv[1];
	stringstream ss;
	ss << argv[2];
	if (!(ss >> radius))
	{
		std::cout << "Error parsing radius(float) " << ss.str() << std::endl;
		exit(-1);
	}
	ss.clear();
	ss << argv[3];
	if (!(ss >> minDist))
	{
		std::cout << "Error parsing minDist(float) " << ss.str() << std::endl;
		exit(-1);
	}
	ss.clear();
	ss << argv[4];
	if (!(ss >> minDens))
	{
		std::cout << "Error parsing minDens(int) " << ss.str() << std::endl;
		exit(-1);
	}
	if (argc >= 6)
	{
		ss.clear();
		ss << argv[5];
		if (!(ss >> std::boolalpha >> exHalo)) {
			std::cout << "Error parsing excludeHalo(bool) " << ss.str() << std::endl;
			exit(-1);
		}
	}
	if (argc >= 7)
	{
		int i;
		ss.clear();
		ss << argv[6];
		if (!(ss >> i))
		{
			std::cout << "Error parsing densType(0/4) " << ss.str() << std::endl;
			exit(-1);
		}
		if (i == 0 || i == 4)
			densType = outputColumn(i);
	}
	if (argc >= 8)
	{
		ss.clear();
		ss << argv[7];
		if (!(ss >> factor))
		{
			std::cout << "Error parsing factor(float >=1) " << ss.str() << std::endl;
			exit(-1);
		}
		if (factor < 1.f)
			factor = 1.f;
	}

	for (int i = 1; i < argc; ++i)
		std::cout << argv[i] << std::endl;

	clustrering(fname, radius, minDist, minDens, exHalo, densType, factor);





	getchar();



	return 0;
}
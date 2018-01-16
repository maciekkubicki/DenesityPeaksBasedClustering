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

#define IF_RANK0 if (world.rank()==0)

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
			//std::cout << "IDX " << i << std::endl;
			output[i][C] = number++;
		}
	return --number;
}

int getClusterNumberFromNeighbour(Matrix<float> output, int neighbour)
{
	cout << neighbour << " ";
	if (output[neighbour][C] != 0)
	{
		cout << "RET " << output[neighbour][C];
		return output[neighbour][C];
	}
	else return getClusterNumberFromNeighbour(output, output[neighbour][N]);
}
std::vector<float> pickClusters(Matrix<float> output, int startIx, int myRows)
{
	std::vector<float> v;
	std::cout << "Picking clusters" << std::endl;

	for (int i = startIx; i < startIx+myRows; ++i)
	{
		std::cout << "pkt " << i << " ";
		if (output[i][C] == 0)
			output[i][C] = getClusterNumberFromNeighbour(output, output[i][N]);
		std::cout << std::endl;
		v.push_back(output[i][C]);
	}
	return v;
}

void normalizeDist(Matrix<float> output, int rows, int cols, float min, float max)
{
	std::cout << "Normalize distances" << std::endl;
	for (int i = 0; i < rows; ++i)
	{
		output[i][D] = (output[i][D] - min) / (max - min);
		output[i][C] = 0;
	}

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
	return sqrt(max);
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
	return sqrt(min);
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
	//std::cout << "Getting min distance from higher denesity" << std::endl;
	float min;
	float max = findMaxDist(mat, rows, cols);
	Matrix<float> q(new float[cols], 1, cols);
	Index<L2<float> > index(mat, flann::KDTreeIndexParams(4));
	index.buildIndex();
	Matrix<int> indices = Matrix<int>(new int[maxR*rows], rows, maxR);
	Matrix<float> dists = Matrix<float>(new float[maxR*rows], rows, maxR);
	index.knnSearch(mat, indices, dists, maxR, flann::SearchParams(128));
	max = maxx(dists, rows, maxR);
	for (int i = 0; i < rows; ++i)
	{
		for (int k = 1; k < maxR; ++k)
		{
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
	std::cout << "Before " << std::endl;
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
	std::cout << "Num of maxes " << _num << "Num of int_max " << _num2 << std::endl;
	if (_num2 < 1)
		return;
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
			q[0][j] = mat[i][j];
		if (output[i][N] == INT_MAX && (i != ix && _num != 1))
		{
			std::cout << "In" << std::endl;
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
				std::cout << "Changed" << std::endl;
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
			if (*rows == 0)
			{
				strline << line;
				while (strline >> test[*cols])
					++*cols;
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
				++j;
			strline.clear();
			++i;
		}
		fp.close();
	}
	std::cout << "RD " << *rows << " " << *cols << std::endl;
	return mat;
}



std::vector<clusterInfo> getMaxBorderDenesity(Matrix<float> output, int numberOfClusters, int rows, outputColumn E = R)
{
	if (E != R && E != CR)
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

void excludeHaloPoints(Matrix<float> output, std::vector<clusterInfo> cInfo, int numberOfClusters, int rows, float factor = 1)
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

void unpackVector(std::vector<std::vector<float>> gatheredVector, Matrix<float> output, outputColumn col, std::vector<int> rowCount, int worldSize)
{
	int c = 0;
	for (int i = 0; i < worldSize; ++i)
		for (int j = 0; j < rowCount[i]; ++j)
			output[c++][col] = gatheredVector[i][j];
}

int takeMyDataset(Matrix<float> dataset, int cols, Matrix<float> myDataset, std::vector<int> rowsPerNode, int myRank, int worldSize)
{
	int startI = 0, endI = 0;
	if (!myRank)
		startI = 0;
	else
		for (int i = 1; i <= myRank; ++i)
			startI += rowsPerNode[i - 1];

	for (int i = 0; i < rowsPerNode[myRank]; ++i)
		for (int j = 0; j < cols; ++j)
			myDataset[i][j] = dataset[i + startI][j];

	return startI;

}
float getMaxiumumDist(Matrix<float> dataset, int rows)
{
	Index<L2<float> > index(dataset, flann::KDTreeIndexParams(4));
	index.buildIndex();
	Matrix<int> indices = Matrix<int>(new int[rows*rows], rows, rows);
	Matrix<float> dists = Matrix<float>(new float[rows*rows], rows, rows);
	index.knnSearch(dataset, indices, dists, rows, flann::SearchParams(128));
	return maxx(dists, rows, rows);
}
void clustrering(boost::mpi::communicator world, const char *fileName, float radius, float minDistF, int minDensF, bool excludeHalo = false, outputColumn densCol = R, float factor = 1.f)
{

	double start = 0, end = 0;
	int numOfClusters = 0;
	int worldSize = world.size();
	int myRank = world.rank();

	std::vector<int> rowsPerNode;


	int a = 2;
	int maxR = 10000;
	int rows, cols;
	//float radius = 4.2;
	float maxDist = 0.f;

	int myRows = 0;
	Matrix<float> dataset;
	IF_RANK0
	{
		dataset = readData(fileName, &rows, &cols);
		rowsPerNode = countRowsPerNode(world.size(), rows);
		maxDist = getMaxiumumDist(dataset, rows);
		start = clock();
		for (int i = 0; i < world.size(); ++i)
			std::cout << "T " << rowsPerNode[i] << "\t";
		std::cout << std::endl;
	}

	std::cout << world.rank() << " afrter reading dataset. Before broadcasting. " << std::endl;
	boost::mpi::broadcast(world, rows, 0);
	boost::mpi::broadcast(world, cols, 0);
	boost::mpi::broadcast(world, dataset, 0);
	boost::mpi::broadcast(world, maxDist, 0);
	boost::mpi::broadcast(world, rowsPerNode, 0);
	boost::mpi::scatter(world, rowsPerNode, myRows, 0);
	world.barrier();

	Matrix<float> extraData(new float[rows * 6], rows, 6);
	Matrix<float> myDataset(new float[myRows*cols], myRows, cols);

	int startIx = takeMyDataset(dataset, cols, myDataset, rowsPerNode, myRank, worldSize);

	Index<L2<float>>index(dataset, flann::LinearIndexParams());
	index.buildIndex();

	IF_RANK0 std::cout << "Start getting denesity" << std::endl;
	
	Matrix<float> q(new float[cols], 1, cols);
	std::vector<float> myExtra;
	std::vector<float> myExtra2(myRows, 0);
	
	for (int i = 0; i < myRows; ++i)
	{
		for (int j = 0; j < cols; ++j)
			q[0][j] = myDataset[i][j];
		Matrix<int> indices = Matrix<int>(new int[maxR], 1, maxR);
		Matrix<float> dists = Matrix<float>(new float[maxR], 1, maxR);
		
		myExtra.push_back(index.radiusSearch(q, indices, dists, radius, flann::SearchParams(128)));
		delete[] indices.ptr();
		delete[] dists.ptr();
	}
	delete[] q.ptr();

	std::cout << "I'm " << world.rank() << " " << rows << " " << cols << " " << myRows << std::endl;
	std::vector<std::vector<float>> all;
	boost::mpi::all_gather(world, myExtra, all);
	world.barrier();
	std::cout << "Tu " << myRank << std::endl;
	IF_RANK0
	{
		unpackVector(all, extraData,R,rowsPerNode,worldSize);
	}
	boost::mpi::broadcast(world, extraData, 0);

	IF_RANK0 std::cout << "Start getting minDist" << std::endl;

	maxR = rows;
	Index<L2<float> > index2(dataset, flann::KDTreeIndexParams(4));
	index2.buildIndex();

	Matrix<int> indices = Matrix<int>(new int[maxR*myRows], myRows, maxR);
	Matrix<float> dists = Matrix<float>(new float[maxR*myRows], myRows, maxR);
	index2.knnSearch(myDataset, indices, dists, maxR, flann::SearchParams(128));
	float max = maxDist;
	float min = std::numeric_limits<float>::min();

	for (int i = 0; i < myRows; ++i)
	{
		if (i % 500 == 0)
			cout << myRank << " " << i << endl;

		for (int k = 1; k < maxR; ++k)
		{
			//cout << i << " " << k << " " << indices[0][k] << " " << dists[0][k] << endl;
			int nn = indices[i][k];
			if (extraData[i + startIx][R] < extraData[nn][R])
			{
				myExtra[i] = dists[i][k];
				myExtra2[i] = nn;
				break;
			}
			else if (k == maxR - 1)
			{

				myExtra[i] = max;
				myExtra2[i] = INT_MAX;
				cout << "Alert " << endl;
				std::cout << i << " " << extraData[i + startIx][R] << " " << k << " " << indices[i][k] << " " << dists[i][k] << endl;

			}
		}

	}
	delete[] indices.ptr();
	delete[] dists.ptr();

	std::vector<std::vector<float>> all2;
	all.clear();
	boost::mpi::all_gather(world, myExtra, all);
	boost::mpi::all_gather(world, myExtra2, all2);
	world.barrier();

	std::cout << "Tu " << myRank << std::endl;
	//float minDistF = 0.02;
	//int minDensF = 250;//rows / 100 * 10;
	IF_RANK0
	{
		unpackVector(all, extraData,D,rowsPerNode,worldSize);
		unpackVector(all2, extraData, N, rowsPerNode, worldSize);
		for (int i = 0; i < rows; ++i)
			std::cout << "all " << i << " " << extraData[i][R] << "\t" << extraData[i][D] << "\t" << extraData[i][N] << std::endl;
		std::cout << std::endl;
		while (countMaxes(extraData, rows) != 1)
		{
			cout << countMaxes(extraData, rows);
			fix2(dataset, extraData, radius, rows, cols, maxR);
		}
		setMinMax(extraData, rows, &min, &max);
		cout << "Got Max " << max << " Min " << min << endl;
		normalizeDist(extraData, rows, cols, min, max);
		numOfClusters = pickClusterCenters(extraData, rows, cols, minDistF, minDensF);
		std::cout << numOfClusters << " Liczba Klastrow" << std::endl;

	}
	boost::mpi::broadcast(world, extraData, 0);
	world.barrier();
	std::vector<float> v = pickClusters(extraData, startIx, myRows);
	all.clear();
	boost::mpi::all_gather(world, v, all);
	world.barrier();
	IF_RANK0 unpackVector(all, extraData,C,rowsPerNode,worldSize);
	boost::mpi::broadcast(world, extraData, 0);
	world.barrier();

	if (excludeHalo)
	{
		q = Matrix<float>(new float[cols], 1, cols);
		//Matrix<float> q(new float[cols], 1, cols);

		for (int i = 0; i < myRows; ++i)
		{
			for (int j = 0; j < cols; ++j)
				q[0][j] = myDataset[i][j];

			Matrix<int> indices = Matrix<int>(new int[maxR], 1, maxR);
			Matrix<float> dists = Matrix<float>(new float[maxR], 1, maxR);

			int nn = index2.radiusSearch(q, indices, dists, radius, flann::SearchParams(128));
			int ix = extraData[startIx + i][C];
			myExtra[i] = 0;
			myExtra2[i] = 0;
			for (int k = 0; k < nn; ++k)
			{
				int ix2 = extraData[indices[0][k]][C];
				if (ix != ix2 && !myExtra2[i])
				{
					std::cout << myRank << " Border point nr " << i + startIx << " bec " << ix << " and " << ix2 << std::endl;
					myExtra[i] = dists[0][k];
					myExtra2[i] = 1;
				}
				else if (ix == ix2)
					++myExtra[i];

			}

		}
		std::cout << myRank << " przed bar" << std::endl;


		all.clear();
		all2.clear();
		boost::mpi::all_gather(world, myExtra, all);
		boost::mpi::all_gather(world, myExtra2, all2);
	}
	world.barrier();

	std::cout << myRank << " po bar" << std::endl;


	IF_RANK0
	{
		if (excludeHalo)
		{
			std::cout << "Unpack" << std::endl;
			unpackVector(all,extraData,CR,rowsPerNode,worldSize);
			std::cout << "Unpack2" << std::endl;
			unpackVector(all2, extraData, iB, rowsPerNode, worldSize);
			std::vector<clusterInfo> ci = getMaxBorderDenesity(extraData, numOfClusters, rows, densCol);
			excludeHaloPoints(extraData, ci, numOfClusters, rows, factor);
		}
		end = clock();
		std::cout << "Pomiar " << std::endl;
		writeData<float>("outputp.dat", extraData, rows, 6);
		double diff = (double)(end - start) / CLOCKS_PER_SEC;
		std::cout << "Time " << diff << std::endl;
	}


	delete[] dataset.ptr();
	delete[] extraData.ptr();
	world.barrier();
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
	boost::mpi::environment env(argc, argv);
	boost::mpi::communicator world;
	IF_RANK0
	{
		std::cout << "Size of world " << world.size() << std::endl;
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
		if(!(ss >> minDens))
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
		
	}
	//boost::mpi::broadcast(world, fname, 0);//only 0
	boost::mpi::broadcast(world, radius, 0);
	boost::mpi::broadcast(world, minDist, 0);
	boost::mpi::broadcast(world, minDens, 0);
	if (argc >= 6) boost::mpi::broadcast(world, exHalo, 0);
	//if (argc >= 7) boost::mpi::broadcast(world, densType, 0);//only 0
	//if (argc >= 8) boost::mpi::broadcast(world, factor, 0);//only 0
	//clustrering(world, "gener2.txt", 4.2, 0.02, 250);
	world.barrier();
	clustrering(world, fname, radius, minDist, minDens, exHalo, densType, factor);
	
	


	return 0;
}


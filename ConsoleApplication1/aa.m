d = importdata('dataset.txt');
d = d(1:4000,:);
dlmwrite('dataset2.txt',d(:,1:4),'delimiter','\t')
k1 = kmeans(d(:,1:4),4);
%d = d(1:100, 2:4);
[coeff, score,latent] = princomp(d);


r = rand(150,4);
%score = score +r;
gscatter(score(:,1), score(:,2), k1(:));
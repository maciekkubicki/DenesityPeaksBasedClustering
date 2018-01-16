clc;
clear all;
 d=iris_dataset';



k1 = kmeans(d(:,1:4),4);

[coeff, score,latent] = princomp(d);


r = rand(150,4);
%score = score +r;
gscatter(score(:,1), score(:,2), k1(:));
%dlmwrite('iris.txt',score(:,1:2),'delimiter','\t','precision',5)

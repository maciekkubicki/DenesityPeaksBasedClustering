d0 = importdata('htru2.txt');
[coeff, score,latent] = princomp(d0);
%score = d0;%%PCA / orginal - comment this line / uncomment this line

k = 2;
[idx,c] = kmeans(d0, k)
d = importdata('outputp.dat');
figure(1)
gscatter(score(:,1), score(:,2), idx);
xlabel('X')
ylabel('Y')
figure(2)
gscatter(score(:,1), score(:,2),d(:,4));
xlabel('X')
ylabel('Y')
figure(3)
histogram(d(:,4))
figure(4)
histogram(idx);

d0 = importdata('gener2.txt');
[coeff, score,latent] = princomp(d0);
score = d0;
k = 7;
[idx,c] = kmeans(d0, k)
d = importdata('outputp.dat');
figure(1)
gscatter(score(:,1), score(:,2), idx);
figure(2)
gscatter(score(:,1), score(:,2),d(:,4));
figure(3)
histogram(d(:,4))
figure(4)
histogram(idx);

d0 = importdata('pcaseeds2.txt');
d = importdata('output.dat');
scatter(d0(:,1), d0(:,2),30,d(:,4));
figure()
histogram(d(:,4))

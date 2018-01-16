o = ones(10,2);
b = randi([-200, 200], 10,2) .* o;
o2 = ones(25,2);
b = [b; randi([20,50], 25,2).*o2];
b = [b; randi([-20,0], 25,2).*o2];
l = randi([150,170],25,1);
p = randi([-50, -20], 25,1);
d = [l,p];
b = [b; d.*o2];

l = randi([10,80],25,1);
p = randi([170, 190], 25,1);
d = [l,p];
b = [b; d.*o2];
scatter(b(:,1), b(:,2));

dlmwrite('myFile.txt',b,'delimiter','\t')

d = importdata('seeds_dataset.txt');
dlmwrite('seeds.txt',d(:,3:4),'delimiter','\t')
figure();
scatter(d(:,3), d(:,4));
x = d(:,3);
y = d(:,4);
max = 0;

for i = 1:size(x,1)
    for j=i+1:size(x,1)
        s = (x(i,1) - x(j,1))^2 + (y(i,1) - y(j,1))^2;
        if(s>max)
            max = s;
        end
    end
end
        

a = 2;
number = 1000;
H = 4; 
N = 1;

V = randi([162, 188], 1, number);
histogram(V)
sum(V)
for i=1:number
    if(V(i)>167 && V(i)<=170)
        V(i)=V(i)-randi([0,5]);
    elseif(V(i)>=170 && V(i)<174)
        V(i) = V(i) + randi([0,5]);
    end
end
figure(2);
sum(V)
histogram(V)
X = [];
Y = [];
for x=162:188
    X = [X, x];
    Y = [Y, estimated(x, H, number, N, V)];
end

figure(3);
plot(X,Y);

[peaks,idx] = findpeaks(Y);
numOfPeaks = size(idx,2);
diffs = [];
for ix = 1:numOfPeaks
    for i = 1:number
        diffs = [diffs, abs(V(i)-162-idx(ix))];
    end
end
XX = [];
YY = [];
for x=0:50
    XX = [XX, x];
    YY = [YY, estimated(x, H, size(diffs,2),N, diffs)];
end

figure(4);
plot(XX,YY);
iYY = max(YY(:)) - YY;
figure(5);
plot(XX,iYY);
[peaks2,idx2] = findpeaks(iYY);
idx2(1)
H = 3;
idx = idx + 162;
vv = zeros(number,2);
for i=1:number
    vv(i, 1) = V(i);
    for ix = 1:numOfPeaks
        if(abs(V(i)-idx(ix))<(idx2(1)-1))
            vv(i,2) = ix;
        end
    end
end
figure(6);
scatter(vv(:,1), vv(:,2));
      
function v = estimated(x,h,m,n,vec)
    v = 0;
    for i=1:m
        v = v + kernel((x-vec(i))/h);
    end
    v = v/(m*h^n);
end

function v = kernel(x)
    if(x>=-1 && x<=1)
        v = 3/4*(1-x^2);
    else 
        v = 0;
    end
end



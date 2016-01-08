function X0 = stand(X)
% STAND将数据矩阵X逐列进行标准化处理，输出标准化数据X0
% X是原始数据矩阵，X0是标准化后的数据矩阵
zeros(size(X));
[nr,nx]=size(X);
for mk=1:nr
X0(mk,:)=(X(mk,:)-mean(X))./std(X);
end
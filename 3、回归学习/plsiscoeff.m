function [COEFF,INTERCEP] = plsiscoeff(X,Y,B)
%PLESTSCOEFF标准化回归系数逆标准化处理，输出原始自变量对因变量的回归系数及常数项
% [COEFF,INTERCEP]=plsiscoeff(X,Y,B)
%X-原始自变量数据
%Y-原始因变量数据
%B-标准化变量回归方程的系数
%COEFF-原始变量回归方程的系数
%INTERCEP-原始变量回归方程的常数项

[xrow,xcol]=size(X);
[yrow,ycol]=size(Y);
for i=1:ycol
    bykCOEFF(:,i)=B(:,i).*std(Y(:,i));
end
for j=1:ycol
    for i=1:xcol
        COEFF(i,j)=bykCOEFF(i,j)/std(X(:,i));
    end
end
INTERCEP=mean(Y)-(mean(X)*COEFF);

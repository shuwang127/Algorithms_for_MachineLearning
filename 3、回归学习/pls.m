function SCOEFF = pls(h,p,W,P,R)
%PLS求偏最小二乘法回归方程的系数
%SCOEFF = pls(h,p,W,P,R)
%h-用于建模的主成分个数
%p-自变量个数
%W-模型效应权重p×rankE0矩阵
%P-模型效应载荷量p×rankE0矩阵
%R-因变量载荷量q×rankE0矩阵
%SCOEFF--偏最小二乘法回归方程的系数p×q矩阵

for byk=1:h
    if byk==1
        WX(:,byk)=W(:,byk);
        SCOEFF=WX(:,byk)*R(:,byk)';
    else
        I=eye(p);
        ww=eye(p);
        for bykbyk=1:byk-1
            ww=ww*(1-W(:,bykbyk)*P(:,bykbyk)');
        end
        WX(:,byk)=ww*W(:,byk);
    end
    SCOEEF=WX(:,byk)*R(:,byk)';
end

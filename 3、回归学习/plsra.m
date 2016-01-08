function RA = plsra(T,R,F0,rankE0)
%PLSRA求出主成分的累积复测定系数
%RA=plsra(T,R,F0,rankE0)
%T—自变量系统主成分得分N×rankE0矩阵
%R—因变量载荷量q×rankE0矩阵
%F0—因变量标准化的样本数据n×q矩阵
%rankE0-plspcr提取的主成分个数

RAAA=[];
for byk=1:rankE0
    RAbyk=sum(norm(T(:,byk)).^2*norm(R(:,byk)).^2)./(norm(F0))^2;
    RAAA=[RAAA,RAbyk];
end
RA=cumsum(RAAA);
    

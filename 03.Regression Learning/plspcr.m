function [W,C,T,U,P,R]=plspcr(E0,F0)
%PLSPCR提取PLS建模过程所有可能的主成分
%[W,C,T,U,P,R]=plspcr(E0,F0)
%E0-自变量标准化的样本数据n×p矩阵；
%F0-因变量标准化的样本数据n×p矩阵；
%W -模拟效应权重p×rankE0矩阵
%C -因变量权重q×rankE0矩阵
%T -自变量系统主成分得分n×rankE0矩阵
%U -因变量系统主成分得分n×rankE0矩阵
%P -模型效应载荷量p×rankE0矩阵
%R -因变量载荷量q×rankE0矩阵

A=rank(E0);
W=[];
C=[];
T=[];
U=[];
P=[];
R=[];
for byk=1:A
%提取主轴与主成分
EFFE=E0'*F0*F0'*E0;
FEEF=F0'*E0*E0'*F0;

[w,LAMBDA]=eigs(EFFE,1,'lm')
[c,LAMBDA]=eigs(FEEF,1,'lm')
t1=E0*w;
u1=F0*c;
W=[W,w];
C=[C,c];
T=[T,t1];
U=[U,u1];

%计算残差
p1=(E0'*t1)/norm(t1)^2;
E1=E0-t1*p1';
E0=E1;
r1=(F0'*t1)/norm(t1)^2;
F1=F0-t1*r1';
F0=F1;
P=[P,p1];
R=[R,r1];
end
















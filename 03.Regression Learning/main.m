%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能：演示回归算法算法在计算机视觉中的应用
%实现如何利用偏最小二乘回归模型实现数据拟合；
%环境：Win7，Matlab2012b
%Modi: NUDT-VAP30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;
%多重相关性计算
load('data.mat'); %读入数据文件data.mat
cr=corrcoef(data) %计算变量之间的相关系数
%建立偏最小二乘回归模型
%提取可能的主成分
X=data(:,1:3);
Y=data(:,3:6);
E0=stand(X)
F0=stand(Y)
A=rank(E0)
[W,C,T,U,P,R]=plspcr(E0,F0)

%主成分解释能力分析
%计算主成分累计复测定系数
RA=plsra(T,R,F0,A)
%计算主成分的信息解释能力
[Rdx,RdX,RdXt,Rdy,RdY,RdYt]=plsrd(E0,F0,T,A)

%计算第一主成分间的相关性
%通过t1/u1图像直观的考查第一主成分间的相关性
cr=plsutcor(U,T)

%计算PLS回归方程的系数
%标准化因变量关于主成分t1的经验回归系数
TCOEFF=R(:,1)
%标准化因变量关于标准化自变量的经验回归系数
SCOEFF=pls(1,5,W,P,R)%1为建模的主成分个数，5为自变量个数
%计算原始因变量关于原始自变量的经验回归系数
[COEFF,INTERCEP]=plsiscoeff(X,Y,SCOEFF)%对标准化的回归系数进行逆标准化处理，输出原始自变量对因变量的回归系数及常数项

%变量投影重要性分析与模型的改进
result=plsresult(W,RdY,RdYt,1)%result表示第j个自变量对因变量的解释能力模
















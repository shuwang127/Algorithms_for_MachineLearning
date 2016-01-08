function [sol] = AcatspEval(sol,distMatrix,numvars)
% 功能－－－计算适应度
% 输入：sol一个个体
% 输出：该个体及其适应度val

val = sum(diag(distMatrix(sol(1:numvars),[sol(2:numvars) sol(1)])));
sol(numvars+1)=val;

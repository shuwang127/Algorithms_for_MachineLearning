%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能：演示随机森林算法在计算机视觉中的应用
%环境：Win7，Matlab2012b
%Modi: NUDT-VAP
%时间：2015-4-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = ['E:\works\book\7（机器学习20讲）\Code\5、Random Forest\'];
data1 = textread([path 'satimage.tra']);
data2 = textread([path 'satimage.txt']);
% path = 'C:\Users\Administrator\Documents\MATLAB\';
% data = textread(path + 'srbct.txt');
% In this data set, each row represents a sample, 
% and each column represents a kind of variable(feature, attribute).
% !! So we should transpose "x" and "xts" below.
[m1, n1] = size(data1);
[m2, n2] = size(data2);  
ntest = m2;  % The number of test set;
ntrain = m1; % The number of training set;
% Above lines we randomly select 2/3 data as training data, 
% and remaining 1/3 data as test data.
x = (data1(1 : ntrain, 1 : n1 - 1));
x = x';
cl = (data1(1 : ntrain, n1));
xts = (data2(1 : ntest, 1 : n2 - 1));
xts = xts';
clts = (data2(1 : ntest, n2));
% Above lines we acquire x, cl, xts and clts from randomData;
nclass = 6;
% The data set has 4 classes.
classwt = 0;
% Here we set all class the same weight 1.
% It can also be written as "classwt = [1 1 1 1];".
cat0 = 0;
% Here we set it having no categorical variables.
runParam = [6 1 50 10 1 0];
% Here we set mtry = 80, ndsize = 1, jbt = 60, look = 10, lookcls = 1, mdim2nd = 0;
impOpt = [0 0 0];
% Here we set imp = 0, Interact = 0, impn = 0;
proCom = [0 0 0 0 0];
% Here we set nprox = 0, nrnn = 0, noutlier = 0, nscale = 0, nprot = 0;
missingVal = 0;
% Here we set missingVal = 0, that means we use the "Default Value" for missingVal.
% That is, code = -999.0, missingfill = 0;
saveForest = [0 0 0];
% Here we set isaverf = 0, isavepar = 0, isavefill = 0;
runForest = [0 0];
% Here we set irunrf = 0, ireadpar = 0;
outParam = [1,0,0,0,0,0,0,0,0,0];
% Here we set isumout = 1 to show a classification summary.
msm = 1 : 36;
% Here we use all 2308 variables, we can also use msm = 0 to use all variables.
seed = 4351;


x = single(x);        %get train x
cl = int32(cl);           %get train label
xts = single(xts);      %get test x
clts = int32(clts);     %get test label
classwt = single(classwt);
cat0 = int32(cat0);
msm = int32(msm);
runParam = int32(runParam);
impOpt = int32(impOpt);
proCom = int32(proCom);
missingVal = single(missingVal);
saveForest = int32(saveForest);
runForest = int32(runForest);
outParam = int32(outParam);
seed = int32(seed);

[errtr, errts, prox, trees, predictts, varimp, scale] = ...
  RF(nclass, x, cl, xts, clts, classwt, cat0, msm, runParam, impOpt, ...
    proCom, missingVal, saveForest, runForest, outParam, seed, 'satimage');
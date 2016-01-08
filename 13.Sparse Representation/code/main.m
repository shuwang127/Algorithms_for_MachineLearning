%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能：演示稀疏表示算法在计算机视觉中的应用
%基于稀疏表示实现人脸识别；
% 实验中采用 SPAMS工具箱和Yale人脸数据库
% SPAMS下载地址 ： 
% http://spams-devel.gforge.inria.fr/
%环境：Win7，Matlab2012b
%Modi: NUDT-VAP
%时间：2014-02-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc, clear all, close all
addpath(genpath('spams-matlab'));
isLoadImg = 0;
%% Load training and testing data
if isLoadImg
dataset = 'Yale';
fileExt = 'png';
datasetPath = ['data/',dataset];
imgDataDir = dir([datasetPath, '*.',fileExt]);
imgSide = 32;

percentTraining = 0.2;
normalization   = true;

fprintf('Loading and sperating training and testing samples ...\n');
[trainLabel, testLabel, trainSample, testSample, numClass] = getTrainTestData(datasetPath, fileExt, imgSide, percentTraining, normalization);

%  Dimension:d --- number of feature 
%            Nt--- number of testing samples
Train.X = trainSample;
Train.y = trainLabel;
Test.X  = testSample;
Test.y  = testLabel;
Nt = size(Test.y,2);
else
    load trainData.mat;
    load testData.mat
    Nt = size(Test.y,2);
end

%% run different methods for classification
% 1. Nearest Neighbor Classifier (NN)
ID1 = knnclassify(Test.X', Train.X', Train.y); % 1-NN Classifier
correctSample = sum(ID1==Test.y');
accuracy(1) = correctSample/Nt;
fprintf('NN - Accuracy = %f %% (%d of %d)\n', accuracy(1)*100, correctSample, Nt);
%----------------------------------------------------------------------------------
% 2. Sparse Representation Classifier (SRC)
% "Face Recognition: A Sparse Representation Perspective" PAMI09
%----------------------------------------------------------------------------------
[ID2, relative_error] = src(Train, Test, 'omp');
correctSample = sum(ID2==Test.y);
accuracy(2) = correctSample/Nt;
fprintf('SRC(OMP Solver) - Accuracy = %f %% (%d of %d)\n', accuracy(2)*100, correctSample, Nt);

[ID3, relative_error] = src(Train, Test, 'apg');
correctSample = sum(ID3==Test.y);
accuracy(3) = correctSample/Nt;
fprintf('SRC(APG Solver) - Accuracy = %f %% (%d of %d)\n', accuracy(3)*100, correctSample, Nt);
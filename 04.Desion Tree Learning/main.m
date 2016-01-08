%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能：演示决策树算法在计算机视觉中的应用
%基于C4.5决策树实现图像二值化；
%环境：Win7，Matlab2012b
%Modi: NUDT-VAP
%时间：2015-4-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all; clear; clc;
%% Step 1  载入图像
Image = imread('flower_test.png');
Mask = imread('flower_mask.png');
figure; imshow(Image); title('Used Image');
figure; imshow(Mask); title('Used Mask');
% In the Mask:
%           Mask(i,j) = 0   -> class 0
%           Mask(i,j) = 255 -> class 1
%           Mask(i,j) = 128 -> unknown
%% Step 2 选择训练数据
[M,N,L] = size(Image);
Data = reshape(Image,[M*N,3]);
pID = find(Mask==255);
nID = find(Mask==0);
pNum = size(pID,1);
nNum = size(nID,1);
% 
TrainData = [Data(pID,:);Data(nID,:)]';
TrainValue = [1*ones([pNum,1]);0*ones([nNum,1])]';
TrainNum = pNum + nNum;
%% Step 3 训练
DivNum = 32;
TrainDataFeatures = uint8(TrainData/DivNum)+1;
Nbins = max(TrainDataFeatures(:));
inc_node = TrainNum*10/100;
discrete_dim = [Nbins,Nbins,Nbins];
tree = BuildC45Tree(TrainDataFeatures, TrainValue, inc_node, discrete_dim, max(discrete_dim));
%% Step 4 测试
TestDataFeatures = uint8(Data'/DivNum)+1;
targets = UseC45Tree(TestDataFeatures, 1:M*N, tree, discrete_dim, unique(TrainValue));
Results = reshape(targets,[M,N]);
% 
figure; imshow(Results,[]); title('C4.5 Classification Results')
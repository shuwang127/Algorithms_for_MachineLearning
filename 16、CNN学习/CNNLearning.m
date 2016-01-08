%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能：演示CNN学习算法在计算机视觉中的应用
%基于CNN网络实现分类；
%环境：Win7，Matlab2012b
%Modi: NUDT-VAP
%时间：2014-10-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load training data
clear;
data = load('../data/mnist.mat');
train_x = double(reshape(data.train_x',28,28,60000))/255;
test_x = double(reshape(data.test_x',28,28,10000))/255;
train_y = double(data.train_y');
test_y = double(data.test_y');
clear('data');

%% CNN setup: set a 5c-2s-10c-2s CNN
cnn.layers = {
    struct('type', 'i') %input layer
    struct('type', 'c', 'outputmaps', 5, 'kernelsize', 5) %convolution layer
    struct('type', 's', 'scale', 2) %sub sampling layer
    struct('type', 'c', 'outputmaps', 10, 'kernelsize', 5) %convolution layer
    struct('type', 's', 'scale', 2) %subsampling layer
};
opts.numepochs      =   2;
opts.batchsize      =  50;
opts.alpha          =   1;
rng('default');
cnn = cnnsetup(cnn, train_x, train_y);

%% CNN train
disp('Start to train CNN:');
[cnn, loss] = cnntrain(cnn, train_x, train_y, opts);

%% CNN test
disp('Start to test CNN:');
[err_rate, err_num] = cnntest(cnn, test_x, test_y);
% With 2 epochs, the error rate: could be around 8%.
disp(['Final classification error rate: ' num2str(err_rate*100), '%.']);

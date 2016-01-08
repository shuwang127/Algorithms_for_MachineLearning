%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能：演示RBM学习算法在计算机视觉中的应用
%RBM学习过程；
%环境：Win7，Matlab2012b
%Modi: NUDT-VAP
%时间：2014-10-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RBMLearning()
%% load training and testing data
clear;
data = load('../data/mnist.mat');
train_x = double(data.train_x) / 255;
clear('data');

%% train RBM model using train_x
% set RBM parameters
num_epochs   = 10;
batch_size   = 100;

input_dim = size(train_x, 2);
hidden_sz = 100;
rbm.alpha    = 1;
rbm.momentum = 0.1;
rbm.W        = zeros(hidden_sz, input_dim);
rbm.vW       = zeros(hidden_sz, input_dim);
rbm.b        = zeros(input_dim, 1);
rbm.vb       = zeros(input_dim, 1);
rbm.c        = zeros(hidden_sz, 1);
rbm.vc       = zeros(hidden_sz, 1);

% train RBM using CD
disp('Start to train RBM:');
rng(0);
m = size(train_x, 1);
num_batches = floor(m / batch_size);
for i = 1 : num_epochs
    kk = randperm(m);
    err = 0;
    for j = 1 : num_batches
        batch = train_x(kk((j-1)*batch_size+1:j*batch_size),:);
        v1 = batch;
        h1 = sigmrnd(repmat(rbm.c', batch_size, 1) + v1 * rbm.W');
        v2 = sigmrnd(repmat(rbm.b', batch_size, 1) + h1 * rbm.W );
        h2 = sigm(repmat(rbm.c', batch_size, 1) + v2 * rbm.W');
        c1 = h1' * v1;
        c2 = h2' * v2;
        
        rbm.vW = rbm.momentum * rbm.vW + rbm.alpha * (c1 - c2) / batch_size;
        rbm.vb = rbm.momentum * rbm.vb + rbm.alpha * sum(v1 - v2)' / batch_size;
        rbm.vc = rbm.momentum * rbm.vc + rbm.alpha * sum(h1 - h2)' / batch_size;
        
        rbm.W = rbm.W + rbm.vW;
        rbm.b = rbm.b + rbm.vb;
        rbm.c = rbm.c + rbm.vc;
        
        err = err + sum(sum((v1 - v2).^ 2)) / batch_size;
    end
    disp(['RBM train: epoch ' num2str(i) '/' num2str(num_epochs)  '. Average reconstruction error is: ' num2str(err / num_batches)]);
end

%% visualize the trained model
V = rbm.W';
minmax = [min(V(:)), max(V(:))];
sz = sqrt(input_dim);
s1 = sz;
s2 = sz;
num=ceil(sqrt(hidden_sz));
data=minmax(2)*ones(num*s2+num-1,num*s1+num-1);
x=0;
y=0;
for i=1:hidden_sz
    im = reshape(V(:,i),s1,s2)';
    data(x*s2+1+x : x*s2+s2+x, y*s1+1+y : y*s1+s1+y)=im;
    x=x+1;
    if(x>=num)
        x=0;
        y=y+1;
    end
end
imagesc(data, [minmax(1) minmax(2)]);
axis equal
axis tight
colormap gray
end
%% other helper functions
function X = sigm(P)
X = 1./(1+exp(-P));
end

function X = sigmrnd(P)
X = double(1./(1+exp(-P)) > rand(size(P)));
end


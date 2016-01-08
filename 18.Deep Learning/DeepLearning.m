%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能：演示深度学习算法在计算机视觉中的应用
%训练DBN用于分类；
%环境：Win7，Matlab2012b
%Modi: NUDT-VAP
%时间：2014-10-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DeepLearning()

%% load training data
clear;
data = load('../data/mnist.mat');
train_x = double(data.train_x) / 255;
train_y = double(data.train_y);
test_x  = double(data.test_x) / 255;
test_y  = double(data.test_y);
clear('data');

%% setup DBN model
input_dim       = size(train_x, 2);
output_dim      = size(train_y, 2);
hidden_sz1      = 100;
hidden_sz2      = 100;
dbn.sizes       = [input_dim, hidden_sz1, hidden_sz2];
dbn_opts.numepochs  =  50;
dbn_opts.batchsize  = 100;
dbn_opts.momentum   = 0.5;
dbn_opts.alpha      =   1;
for i = 1:numel(dbn.sizes)-1
    dbn.rbm{i}.alpha    = dbn_opts.alpha;
    dbn.rbm{i}.momentum = dbn_opts.momentum;
    dbn.rbm{i}.W        = zeros(dbn.sizes(i+1), dbn.sizes(i));
    dbn.rbm{i}.vW       = zeros(dbn.sizes(i+1), dbn.sizes(i));
    dbn.rbm{i}.b        = zeros(dbn.sizes(i), 1);
    dbn.rbm{i}.vb       = zeros(dbn.sizes(i), 1);
    dbn.rbm{i}.c        = zeros(dbn.sizes(i+1), 1);
    dbn.rbm{i}.vc       = zeros(dbn.sizes(i+1), 1);
end

%% layer-wise pretrain DBN model
x = train_x;
rng('default');
dbn.rbm{1} = rbmtrain(dbn.rbm{1}, x, dbn_opts);
for i = 2:numel(dbn.sizes)-1
    x = sigm(repmat(dbn.rbm{i-1}.c', size(x, 1), 1) + x * dbn.rbm{i-1}.W');
    dbn.rbm{i} = rbmtrain(dbn.rbm{i}, x, dbn_opts);
end

%% use DBN to initialize NN to peform fine-tuning
nn.sizes            = [dbn.sizes, output_dim];
nn.n                = numel(nn.sizes);
nn.learning_rate    = 0.1;
nn.momentum         = 0.5;

nn_opts.numepochs   =  50;
nn_opts.batchsize   = 100;


rng('default');
for i = 1:numel(dbn.rbm)
    nn.W{i}     = [dbn.rbm{i}.c, dbn.rbm{i}.W];
    nn.vW{i}    = zeros(size(nn.W{i}));
end
i = nn.n - 1;
nn.W{i}  = (rand(nn.sizes(i+1), nn.sizes(i)+1) - 0.5) * 2 * 4 * sqrt(6 / (nn.sizes(i+1) + nn.sizes(i)));
nn.vW{i} = zeros(size(nn.W{i}));

nn = nntrain(nn, train_x, train_y, nn_opts);
[err_rate, ~] = nntest(nn, test_x, test_y);

% With 50 epochs, the error rate: could be around 5%.
disp(['Final classification error rate: ' num2str(err_rate*100), '%.']);

end

function X = sigm(P)
    X = 1./(1+exp(-P));
end












function [net, L] = cnntrain(net, x, y, opts)
m = size(x, 3);
numbatches = m / opts.batchsize;
if rem(numbatches, 1) ~= 0
    error('numbatches not integer');
end
L = zeros(opts.numepochs*numbatches,1);
n  = 1;
for i = 1 : opts.numepochs
    tic;
    kk = randperm(m);
    for k = 1 : numbatches
        batch_x = x(:, :, kk((k - 1) * opts.batchsize + 1 : k * opts.batchsize));
        batch_y = y(:,    kk((k - 1) * opts.batchsize + 1 : k * opts.batchsize));
        
        net = cnnff(net, batch_x);
        net = cnnbp(net, batch_y);
        net = cnngrads(net, opts);
        
        L(n) = net.loss;
        n = n + 1;
    end
    t = toc;
    str_perf = sprintf('; Full-batch train err = %f', net.loss);
    disp(['CNN train: epoch ' num2str(i) '/' num2str(opts.numepochs) '. Took ' num2str(t) ' seconds' '. Mini-batch mean squared error on training set is ' num2str(mean(L((n-numbatches):(n-1)))) str_perf]);
end
end

function net = cnnff(net, x)
n = numel(net.layers);
net.layers{1}.a{1} = x;
inputmaps = 1;

for l = 2 : n   %  for each layer
    if strcmp(net.layers{l}.type, 'c')
        %  !!below can probably be handled by insane matrix operations
        for j = 1 : net.layers{l}.outputmaps   %  for each output map
            %  create temp output map
            z = zeros(size(net.layers{l - 1}.a{1}) - [net.layers{l}.kernelsize - 1 net.layers{l}.kernelsize - 1 0]);
            for i = 1 : inputmaps   %  for each input map
                %  convolve with corresponding kernel and add to temp output map
                z = z + convn(net.layers{l - 1}.a{i}, net.layers{l}.k{i}{j}, 'valid');
            end
            %  add bias, pass through nonlinearity
            net.layers{l}.a{j} = sigm(z + net.layers{l}.b{j});
        end
        %  set number of input maps to this layers number of outputmaps
        inputmaps = net.layers{l}.outputmaps;
    elseif strcmp(net.layers{l}.type, 's')
        %  downsample
        for j = 1 : inputmaps
            z = convn(net.layers{l - 1}.a{j}, ones(net.layers{l}.scale) / (net.layers{l}.scale ^ 2), 'valid');   %  !! replace with variable
            net.layers{l}.a{j} = z(1 : net.layers{l}.scale : end, 1 : net.layers{l}.scale : end, :);
        end
    end
end

%  concatenate all end layer feature maps into vector
net.fv = [];
for j = 1 : numel(net.layers{n}.a)
    sa = size(net.layers{n}.a{j});
    net.fv = [net.fv; reshape(net.layers{n}.a{j}, sa(1) * sa(2), sa(3))];
end
%  feedforward into output perceptrons
net.o = sigm(net.ffW * net.fv + repmat(net.ffb, 1, size(net.fv, 2)));
end

function X = sigm(P)
X = 1./(1+exp(-P));
end

function net = cnnbp(net, y)
n = numel(net.layers);

%   error
net.e = net.o - y;
%  loss function
net.loss = 1/2* sum(net.e(:) .^ 2) / size(net.e, 2);

%%  backprop deltas
net.od = net.e .* (net.o .* (1 - net.o));   %  output delta
net.fvd = (net.ffW' * net.od);              %  feature vector delta
if strcmp(net.layers{n}.type, 'c')         %  only conv layers has sigm function
    net.fvd = net.fvd .* (net.fv .* (1 - net.fv));
end

%  reshape feature vector deltas into output map style
sa = size(net.layers{n}.a{1});
fvnum = sa(1) * sa(2);
for j = 1 : numel(net.layers{n}.a)
    net.layers{n}.d{j} = reshape(net.fvd(((j - 1) * fvnum + 1) : j * fvnum, :), sa(1), sa(2), sa(3));
end

for l = (n - 1) : -1 : 1
    if strcmp(net.layers{l}.type, 'c')
        for j = 1 : numel(net.layers{l}.a)
            net.layers{l}.d{j} = net.layers{l}.a{j} .* (1 - net.layers{l}.a{j}) .* (expand(net.layers{l + 1}.d{j}, [net.layers{l + 1}.scale net.layers{l + 1}.scale 1]) / net.layers{l + 1}.scale ^ 2);
        end
    elseif strcmp(net.layers{l}.type, 's')
        for i = 1 : numel(net.layers{l}.a)
            z = zeros(size(net.layers{l}.a{1}));
            for j = 1 : numel(net.layers{l + 1}.a)
                z = z + convn(net.layers{l + 1}.d{j}, rot180(net.layers{l + 1}.k{i}{j}), 'full');
            end
            net.layers{l}.d{i} = z;
        end
    end
end

%%  calc gradients
for l = 2 : n
    if strcmp(net.layers{l}.type, 'c')
        for j = 1 : numel(net.layers{l}.a)
            for i = 1 : numel(net.layers{l - 1}.a)
                net.layers{l}.dk{i}{j} = convn(flipall(net.layers{l - 1}.a{i}), net.layers{l}.d{j}, 'valid') / size(net.layers{l}.d{j}, 3);
            end
            net.layers{l}.db{j} = sum(net.layers{l}.d{j}(:)) / size(net.layers{l}.d{j}, 3);
        end
    end
end
net.dffW = net.od * (net.fv)' / size(net.od, 2);
net.dffb = mean(net.od, 2);
end

function X = rot180(X)
X = flip(flip(X, 1), 2);
end

function B = expand(A, S)
SA = size(A);  % Get the size (and number of dimensions) of input.
if length(SA) ~= length(S)
    error('Length of size vector must equal ndims(A).  See help.')
elseif any(S ~= floor(S))
    error('The size vector must contain integers only.  See help.')
end
T = cell(length(SA), 1);
for ii = length(SA) : -1 : 1
    H = zeros(SA(ii) * S(ii), 1);   %  One index vector into A for each dim.
    H(1 : S(ii) : SA(ii) * S(ii)) = 1;   %  Put ones in correct places.
    T{ii} = cumsum(H);   %  Cumsumming creates the correct order.
end
B = A(T{:});   %  Feed the indices into A.
end

function X=flipall(X)
for i=1:ndims(X)
    X = flip(X,i);
end
end



function net = cnngrads(net, opts)
for l = 2 : numel(net.layers)
    if strcmp(net.layers{l}.type, 'c')
        for j = 1 : numel(net.layers{l}.a)
            for ii = 1 : numel(net.layers{l - 1}.a)
                net.layers{l}.k{ii}{j} = net.layers{l}.k{ii}{j} - opts.alpha * net.layers{l}.dk{ii}{j};
            end
            net.layers{l}.b{j} = net.layers{l}.b{j} - opts.alpha * net.layers{l}.db{j};
        end
    end
end

net.ffW = net.ffW - opts.alpha * net.dffW;
net.ffb = net.ffb - opts.alpha * net.dffb;
end

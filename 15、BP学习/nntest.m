function [err_rate, err_num] = nntest(nn, x, y)
nn = nnff(nn, x, zeros(size(x,1), nn.sizes(end)));
[~, labels] = max(nn.a{end},[],2);
[~, expected] = max(y,[],2);
err_num = find(labels ~= expected);
err_rate = numel(err_num) / size(x, 1);
end

function nn = nnff(nn, x, y)
n = nn.n;
m = size(x, 1);
x = [ones(m,1) x];
nn.a{1} = x;
%feedforward pass
for i = 2 : n-1
    nn.a{i} = sigm(nn.a{i - 1} * nn.W{i - 1}');
    %Add the bias term
    nn.a{i} = [ones(m,1) nn.a{i}];
end
nn.a{n} = sigm(nn.a{n - 1} * nn.W{n - 1}');
%error and loss
nn.error = y - nn.a{n};
nn.loss = 1/2 * sum(sum(nn.error .^ 2)) / m;

function X = sigm(P)
    X = 1./(1+exp(-P));
end

end

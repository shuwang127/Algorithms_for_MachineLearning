function [K] = kernel(ker,x,y)
% Calculate kernel function.   

switch ker.type
    case 'linear'
        K = x'*y;
    case 'ploy'
        d = ker.degree;
        c = ker.offset;
        K = (x'*y+c).^d;
    case 'gauss'
        
        s = ker.width;
        rows = size(x,2);
        cols = size(y,2);   
        tmp = zeros(rows,cols);
        for i = 1:rows
            for j = 1:cols
                tmp(i,j) = norm(x(:,i)-y(:,j));
            end
        end        
        K = exp(-0.5*(tmp/s).^2);

    case 'tanh'
        g = ker.gamma;
        c = ker.offset;
        K = tanh(g*x'*y+c);
    otherwise
        K = 0;
end
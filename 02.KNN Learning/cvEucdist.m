% cvEucdist - Euclidean distance
%
% Synopsis
%   [d] = cvEucdist(X, Y)
%
% Description
%   cvEucdist calculates a squared euclidean distance between X and Y.
%
% Inputs ([]s are optional)
%   (matrix) X        D x N matrix where D is the dimension of vectors
%                     and N is the number of vectors.
%   (matrix) [Y]      D x P matrix where D is the dimension of vectors
%                     and P is the number of vectors.
%                     If Y is not given, the L2 norm of X is computed and
%                     1 x N matrix (not N x 1) is returned.
%
% Outputs ([]s are optional)
%   (matrix) d        N x P matrix where d(n,p) represents the squared
%                     euclidean distance between X(:,n) and Y(:,p).
%
% Examples
%   X = [1 2
%        1 2];
%   Y = [1 2 3
%        1 2 3];
%   d = cvEucdist(X, Y)
% %      0     2     8
% %      2     0     2
%
% See also
%   cvMahaldist

% Authors
%   Naotoshi Seo <sonots(at)sonots.com>
%
% License
%   The program is free to use for non-commercial academic purposes,
%   but for course works, you must understand what is going inside to use.
%   The program can be used, modified, or re-distributed for any purposes
%   if you or one of your group understand codes (the one must come to
%   court if court cases occur.) Please contact the authors if you are
%   interested in using the program without meeting the above conditions.
%
% Changes
%   06/2006  First Edition
function d = cvEucdist(X, Y)
 if ~exist('Y', 'var') || isempty(Y)
     %% Y = zeros(size(X, 1), 1);
     U = ones(size(X, 1), 1);
     d = abs(X'.^2*U).'; return;
 end
 V = ~isnan(X); X(~V) = 0; % V = ones(D, N); 
 %clear V;
 U = ~isnan(Y); Y(~U) = 0; % U = ones(D, P); 
 %clear U;
 %d = abs(X'.^2*U - 2*X'*Y + V'*Y.^2);
 d1 = X'.^2*U;
 d3 = V'*Y.^2;
 d2 = X'*Y;
 d = abs(d1-2*d2+d3);
 
% X = X';
% Y = Y';
% for i=1:size(X,1)
%     for j=1:size(Y,1)
%         d(i,j)=(norm(X(i,:)-Y(j,:)))^2;  %计算每个测试样本与所有训练样本的欧氏距离
%     end
% end




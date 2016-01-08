% cvKnn - K-Nearest Neighbor classification
%
% Synopsis
%   [Class] = cvKnn(X, Proto, ProtoClass, [K], [distFunc])
%
% Description
%   K-Nearest Neighbor classification
%
% Inputs ([]s are optional)
%   (matrix) X        D x N matrix representing column classifiee vectors
%                     where D is the number of dimensions and N is the
%                     number of vectors.
%   (matrix) Proto    D x P matrix representing column prototype vectors
%                     where D is the number of dimensions and P is the
%                     number of vectors.
%   (vector) ProtoClass
%                     1 x P vector containing class lables for prototype
%                     vectors. 
%   (scalar) [K = 1]  K-NN's K. Search K nearest neighbors.
%   (func)   [distFunc = @cvEucdist]
%                     A function handle for distance measure. The function
%                     must have two arguments for matrix X and Y. See
%                     cvEucdist.m (Euclidean distance) as a reference.
%
% Outputs ([]s are optional)
%   (vector) Class    1 x N vector containing classified class labels 
%                     for X. Class(n) is the class id for X(:,n). 
%   (matrix) [Rank]   Available only for NN (K = 1) now.
%                     nClass x N vector containing ranking class labels
%                     for X. Rank(1,n) is the 1st candidate which is 
%                     the same with Class(n), Rank(2,n) is the 2nd 
%                     candidate, Rank(3,n) is the 3rd, and so on.
%
% See also
%   cvEucdist, cvMahaldist

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
%   04/01/2005  First Edition
function [Class, Rank] = cvKnn(X, Proto, ProtoClass, K, distFunc)
if ~exist('K', 'var') || isempty(K)
    K = 1;%默认为K = 1
end
if ~exist('distFunc', 'var') || isempty(distFunc)
    distFunc = @cvEucdist;
end
if size(X, 1) ~= size(Proto, 1)
    error('Dimensions of classifiee vectors and prototype vectors do not match.');
end
[D, N] = size(X);

% Calculate euclidean distances between classifiees and prototypes
d = distFunc(X, Proto);

if K == 1, % sort distances only if K>1
    [mini, IndexProto] = min(d, [], 2); % 2 == row%每列的最小元素
    Class = ProtoClass(IndexProto);
    if nargout == 2, % instance indices in similarity descending order
        [sorted, ind] = sort(d'); % PxN
        RankIndex = ProtoClass(ind); %,e.g., [2 1 2 3 1 5 4 1 2]'
        % conv into, e.g., [2 1 3 5 4]'
        for n = 1:N
            [ClassLabel, ind] = unique(RankIndex(:,n),'first');
            [sorted, ind] = sort(ind);
            Rank(:,n) = ClassLabel(ind);
        end
    end
else
    [sorted, IndexProto] = sort(d'); % PxN
    clear d；
    % K closest
    IndexProto = IndexProto(1:K,:);
    KnnClass = ProtoClass(IndexProto);
    % Find all class labels
    ClassLabel = unique(ProtoClass);
    nClass = length(ClassLabel);
    for i = 1:nClass
        ClassCounter(i,:) = sum(KnnClass == ClassLabel(i));
    end
    [maxi, winnerLabelIndex] = max(ClassCounter, [], 1); % 1 == col
    % Future Work: Handle ties somehow
    Class = ClassLabel(winnerLabelIndex);
end
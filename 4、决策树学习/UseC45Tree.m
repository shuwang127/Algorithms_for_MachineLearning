function targets = UseC45Tree(features, indices, tree, discrete_dim, Uc)
% Input:
%       features:       Ni*L, Ni features and L samples
%       indices:        index
%       tree:           C4.5 decision tree
%       discrete_dim:   discrete number in each dim
%       Uc:             target classes
% Output:
%       targets:        classification results
%% Step 0-initilize the results
targets = zeros(1, size(features,2));
%% Step 1-Stop condition: leaf node
if (tree.dim == 0)
   % Reached the end of the tree
   targets(indices) = tree.child;
   return;
end
%% Step 2-Otherwise: use children node
% 2-1 First, find the dimension we are to work on
dim = tree.dim;
dims= 1:size(features,1);
% 2-2 And classify according to it
if (discrete_dim(dim) == 0),
    % Continuous feature
    in			= indices(find(features(dim, indices) <= tree.split_loc));
    targets		= targets + UseC45Tree(features(dims, :), in, tree.child(1), discrete_dim(dims), Uc);
    in			= indices(find(features(dim, indices) >  tree.split_loc));
    targets		= targets + UseC45Tree(features(dims, :), in, tree.child(2), discrete_dim(dims), Uc);
else
    % Discrete feature
    Uf			= unique(features(dim,:));
    for i = 1:length(Uf),
        in   	= indices(find(features(dim, indices) == Uf(i)));
        targets	= targets + UseC45Tree(features(dims, :), in, tree.child(i), discrete_dim(dims), Uc);
   end
end
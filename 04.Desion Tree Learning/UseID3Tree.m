function targets = UseID3Tree(features, indices, tree, Nbins, Uc)
% Input:
%       features:   Ni*L, Ni features and L samples
%       indices:    1:Num
%       tree:       ID3 tree
%       Nbins:      features are in 1:Nbins
%       Uc:         target class ID
% Output:
%       targets:    classification results
%% Step 0-initilize the results
targets = zeros(1, size(features,2));
%% Step 1-Stop Condition: if the feature is one dimension
if (size(features,1) == 1)
    for i = 1:Nbins
        in = indices(find(features(indices) == i));
        if ~isempty(in)
            if isfinite(tree.child(i))
                targets(in) = tree.child(i);
            else
                % No data was found in the training set for this bin, so choose it randomally
                n           = 1 + floor(rand(1)*length(Uc));
                targets(in) = Uc(n);
            end
        end
    end
    return;
end
%% Step 2-Otherwise: use the children node for classification
% 2-1 First, find the dimension we are to work on
dim = tree.split_dim;
dims= find(~ismember(1:size(features,1), dim));
% 2-2 And classify according to it
for i = 1:Nbins
    in      = indices(find(features(dim, indices) == i));
    targets = targets + UseID3Tree(features(dims, :), in, tree.child(i), Nbins, Uc);
end
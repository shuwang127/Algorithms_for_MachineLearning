function tree = BuildC45Tree(features, targets, inc_node, discrete_dim, maxNbin)
% Input:
%       features:       Ni*L, Ni features and L samples
%       targets:        Uc classes
%       inc_node:       allow incorrect number 
%       discrete_dim:   discrete number in each dim
%       maxNbin:        max(discrete_dim)
%       base:           default = 0
% Output:
%       tree:           C4.5 decision tree
%% Step 0-Get the size
[Ni, L]         = size(features);
Uc              = unique(targets);
% set default value
tree.dim        = 0;
% tree.child(1:maxNbin)	= zeros(1,maxNbin);
tree.split_loc  = inf;
if isempty(features),
    return;
end
%% Step 1-Stop Condition: feature dim is one or examples is small
if ((inc_node > L) || (L == 1) || (length(Uc) == 1)),
    H = hist(targets, length(Uc));
    [m, largest] 	= max(H);
    tree.child	 	= Uc(largest);
    return;
end
%% Step 2-Otherwise: use C4.5 choose the best feature
% 2-1 Compute the node's I
for i = 1:length(Uc),
    Pnode(i) = length(find(targets == Uc(i))) / L;
end
Inode = -sum(Pnode.*log(Pnode)/log(2));
% 2-2 For each dimension, compute the gain ratio impurity
% This is done separately for discrete and continuous features
delta_Ib    = zeros(1, Ni);
split_loc	= ones(1, Ni)*inf;
for i = 1:Ni,
    data	= features(i,:);
    Nbins	= length(unique(data));
    if (discrete_dim(i)),
        %This is a discrete feature
        P	= zeros(length(Uc), Nbins);
        for j = 1:length(Uc),
            for k = 1:Nbins,
                indices = find((targets == Uc(j)) & (features(i,:) == k));
                P(j,k)  = length(indices);
            end
        end
        Pk	= sum(P);
        P   = P/L;
        Pk          = Pk/sum(Pk);
        info        = sum(-P.*log(eps+P)/log(2));
        delta_Ib(i) = (Inode-sum(Pk.*info))/-sum(Pk.*log(eps+Pk)/log(2));
    else
        % This is a continuous feature
        P	= zeros(length(Uc), 2);
        % Sort the features
        [sorted_data, indices] = sort(data);
        sorted_targets = targets(indices);
        % Calculate the information for each possible split
        I	= zeros(1, L-1);
        for j = 1:L-1,
            for k =1:length(Uc),
                P(k,1) = length(find(sorted_targets(1:j) 		== Uc(k)));
                P(k,2) = length(find(sorted_targets(j+1:end) == Uc(k)));
            end
            Ps		= sum(P)/L;
            P		= P/L;
            info	= sum(-P.*log(eps+P)/log(2));
            I(j)	= Inode - sum(info.*Ps);   
        end
        [delta_Ib(i), s] = max(I);
		split_loc(i) = sorted_data(s);      
    end
end
% 2-3 Find the dimension minimizing delta_Ib 
[m, dim] = max(delta_Ib);
dims     = 1:Ni;
tree.dim = dim;
% 2-4 Split along the 'dim' dimension
Nf       = unique(features(dim,:));
Nbins	 = length(Nf);
if (discrete_dim(dim)),
    %Discrete feature
    for i = 1:Nbins,
        indices    		= find(features(dim, :) == Nf(i));
        tree.child(i)	= BuildC45Tree(features(dims, indices), targets(indices), inc_node, discrete_dim(dims), maxNbin);
    end
else
    %Continuous feature
    tree.split_loc		= split_loc(dim);
    indices1		   	= find(features(dim,:) <= split_loc(dim));
    indices2	   		= find(features(dim,:) > split_loc(dim));
    tree.child(1)		= BuildC45Tree(features(dims, indices1), targets(indices1), inc_node, discrete_dim(dims), maxNbin);
    tree.child(2)		= BuildC45Tree(features(dims, indices2), targets(indices2), inc_node, discrete_dim(dims), maxNbin);
end
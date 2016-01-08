%   This file is part of GML Matlab Toolbox
%   For conditions of distribution and use, see the accompanying License.txt file.
%
%   do_learn_nu Implements splitting of tree node
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%    [tree_node_left, tree_node_right, split_error] =
%    do_learn_nu(tree_node, dataset, labels, weights, papa)
%    ---------------------------------------------------------------------------------
%    Arguments:
%           tree_node - object of tree_node_w class
%           dataset   - training data
%           labels    - training labels
%           weights   - weights of training data
%           papa      - parent node (the one being split)
%    Return:
%           tree_node_left  - left node (result of splitting)
%           tree_node_right - right node (result of splitting)
%           split_error     - error of splitting

function [tree_node_left, tree_node_right, split_error] = do_learn_nu(tree_node, dataset, labels, weights, papa)

tree_node_left = tree_node;
tree_node_right = tree_node;

if(nargin > 4)
  tree_node_left.parent  = papa;
  tree_node_right.parent = papa;
end

Distr = weights;

%[trainpat, traintarg] = get_train( dataset);
trainpat = dataset;
traintarg = labels;

tr_size = size(trainpat, 2);

T_MIN = zeros(3,size(trainpat,1));
d_min = 1;
d_max = size(trainpat,1);

for d = d_min : d_max;

  [DS, IX] = sort(trainpat(d,:));

  TS = traintarg(IX);
  DiS = Distr(IX);
    
  lDS = length(DS);
  
  vPos = 0 * TS;
  vNeg = vPos;
  
  i = 1;
  j = 1;
  
  while i <= lDS
    k = 0;
    while i + k <= lDS && DS(i) == DS(i+k)
      if(TS(i+k) > 0)
        vPos(j) = vPos(j) + DiS(i+k);
      else
        vNeg(j) = vNeg(j) + DiS(i+k);
      end
      k = k + 1;
    end
    i = i + k;
    j = j + 1;
  end
  
  vNeg = vNeg(1:j-1);
  vPos = vPos(1:j-1);
  
  Error = zeros(1, j - 1);

  InvError = Error;
  
  IPos = vPos;
  INeg = vNeg;
  
  for i = 2 : length(IPos)
    IPos(i) = IPos(i-1) + vPos(i);
    INeg(i) = INeg(i-1) + vNeg(i);
  end
  
  Ntot = INeg(end);
  Ptot = IPos(end);
  
  for i = 1 : j - 1
    Error(i) = IPos(i) + Ntot - INeg(i);
    InvError(i) = INeg(i) + Ptot - IPos(i);
  end
  
  idx_of_err_min = find(Error == min(Error));
  if(length(idx_of_err_min) < 1)
      idx_of_err_min = 1;  
  end
  
  if(length(idx_of_err_min) <1)
    idx_of_err_min = idx_of_err_min;
  end
  idx_of_err_min = idx_of_err_min(1);
  
  idx_of_inv_err_min = find(InvError == min(InvError));
  
  if(length(idx_of_inv_err_min) < 1)
      idx_of_inv_err_min = 1;  
  end
  
  idx_of_inv_err_min = idx_of_inv_err_min(1);
  
  if(Error(idx_of_err_min) < InvError(idx_of_inv_err_min))
    T_MIN(1,d) = Error(idx_of_err_min);
    T_MIN(2,d) = idx_of_err_min;
    T_MIN(3,d) = -1;
  else
    T_MIN(1,d) = InvError(idx_of_inv_err_min);
    T_MIN(2,d) = idx_of_inv_err_min;
    T_MIN(3,d) = 1;
  end
  
end

dim = [];

best_dim = find(T_MIN(1,:) == min(T_MIN(1,:)));

dim = best_dim(1);  

tree_node_left.dim = dim;
tree_node_right.dim = dim;

TDS = sort(trainpat(dim,:));

lDS = length(TDS);

DS = TDS * 0;

i = 1;
j = 1;

while i <= lDS
  k = 0;
  while i + k <= lDS && TDS(i) == TDS(i+k) 
    DS(j) = TDS(i);
    k = k + 1;
  end
  i = i + k;
  j = j + 1;
end

DS = DS(1:j-1);

split = (DS(T_MIN(2,dim)) + DS(min(T_MIN(2,dim) + 1, length(DS)))) / 2;

split_error = T_MIN(1,dim);

tree_node_left.right_constrain = split;
tree_node_right.left_constrain = split;

function [i,t] = weakLearner(distribution,train,label)
	%disp('run weakLearner');
    for tt = unique(train)%1:(16*256-1)    
        error(tt)=(label .* distribution) * ((train(:,floor(tt/16)+1)>=16*(mod(tt,16)+1)));         
    end
    [val,tt]=max(abs(error-0.5));
    
    i=floor(tt/16)+1;
    t=16*(mod(tt,16)+1);
    return;
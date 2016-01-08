%   This file is part of GML Matlab Toolbox
%   For conditions of distribution and use, see the accompanying License.txt file.
%
%   train Implements training of a classification tree
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%    nodes = train(node, dataset, labels, weights)
%    ---------------------------------------------------------------------------------
%    Arguments:
%           node      - object of tree_node_w class (initialized properly)
%           dataset   - training data
%           labels    - training labels
%           weights   - weights of training data
%    Return:
%           nodes     - tree is represented as a cell array of its nodes

function nodes = train(node, dataset, labels, weights)

max_split = node.max_split;

[left right spit_error] = do_learn_nu(node, dataset, labels, weights);

nodes = {left, right};

left_pos  = sum((calc_output(left , dataset) == labels) .* weights);
left_neg  = sum((calc_output(left , dataset) == -labels) .* weights);
right_pos = sum((calc_output(right, dataset) == labels) .* weights);
right_neg = sum((calc_output(right, dataset) == -labels) .* weights);

errors = [min(left_pos, left_neg), min(right_pos, right_neg)];

if(right_pos == 0 && right_neg == 0)
  return;
end

if(left_pos == 0 && left_neg == 0)
  return;
end

[errors, IDX] = sort(errors);
errors = flipdim(errors,2);
IDX    = flipdim(IDX,2);
nodes  = nodes(IDX);


splits = [];
split_errors = [];
deltas = [];


for i = 2 : max_split
    for j = 1 : length(errors)
        
        if(length(deltas) >= j)
            continue;
        end
        
        max_node = nodes{j};
        max_node_out = calc_output(max_node, dataset);
       
        mask = find(max_node_out == 1);  
       
        [left right spit_error] = do_learn_nu(node, dataset(:,mask), labels(mask), weights(mask), max_node);
              
        
        left_pos  = sum((calc_output(left , dataset) == labels) .* weights);
        left_neg  = sum((calc_output(left , dataset) == -labels) .* weights);
        right_pos = sum((calc_output(right, dataset) == labels) .* weights);
        right_neg = sum((calc_output(right, dataset) == -labels) .* weights);
        
        splits{end+1} = left;
        splits{end+1} = right;  
        
        if( (right_pos + right_neg) == 0 || (left_pos + left_neg) == 0)
          deltas(end+1) = 0;
        else
          deltas(end+1) = errors(j) - spit_error;
        end
        
        split_errors(end+1) = min(left_pos, left_neg);
        split_errors(end+1) = min(right_pos, right_neg);
    end  
    
    if(max(deltas) == 0)
        return;
    end
    best_split = find(deltas == max(deltas));
    best_split = best_split(1);
    
    cut_vec = [1 : (best_split-1)  (best_split + 1) : length(errors)];
    nodes   = nodes(cut_vec);
    errors  = errors(cut_vec);
    deltas  = deltas(cut_vec);
    
    nodes{end+1} = splits{2 * best_split - 1};
    nodes{end+1} = splits{2 * best_split};
    
    errors(end+1) = split_errors(2 * best_split - 1);
    errors(end+1) = split_errors(2 * best_split);
    
    cut_vec = [1 : 2 * (best_split-1)  2 * (best_split)+1 : length(split_errors)];
    split_errors = split_errors(cut_vec);    
    splits       = splits(cut_vec);

end
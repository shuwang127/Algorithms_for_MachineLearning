%   This file is part of GML Matlab Toolbox
%   For conditions of distribution and use, see the accompanying License.txt file.
%
%   get_dim_and_tr is the function, that returns dimension and threshold of
%   tree node
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%    output = get_dim_and_tr(tree_node, output)
%    ---------------------------------------------------------------------------------
%    Arguments:
%           tree_node - a node of classification tree
%           output    - vector of dimensions and thresholds. Result fo
%                       current node would be concatinated to it
%    Return:
%           output    - a vector of thresholds and dimensions. It has the
%                       following format:
%                       [dimension threshold left/right ...]
%                       left/right is  [-1, +1] number, wich signifies if
%                       current threshold is eather left or right

function output = get_dim_and_tr(tree_node, output)

if(nargin < 2)
  output = [];
end

if(length(tree_node.parent) > 0)
  output = get_dim_and_tr(tree_node.parent, output);
end

output(end+1) = tree_node.dim;

if( length(tree_node.right_constrain) > 0)
  output(end+1) = tree_node.right_constrain;
  output(end+1) = -1;
elseif( length(tree_node.left_constrain) > 0)
  output(end+1) = tree_node.left_constrain;
  output(end+1) = +1;
end


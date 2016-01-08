%   This file is part of GML Matlab Toolbox
%   For conditions of distribution and use, see the accompanying License.txt file.
%
%   calc_output Implements classification of input by a classification tree node
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%    y = calc_output(tree_node, XData)
%    ---------------------------------------------------------------------------------
%    Arguments:
%           tree_node - classification tree node
%           XData     - data, that will be classified
%    Return:
%           y         - +1, if XData belongs to tree node, -1 otherwise (y is a vector)

function y = calc_output(tree_node, XData)
y = XData(tree_node.dim, :) * 0 + 1;


for i = 1 : length(tree_node.parent)
  y = y .* calc_output(tree_node.parent, XData);
end

if( length(tree_node.right_constrain) > 0)
  y = y .* ((XData(tree_node.dim, :) < tree_node.right_constrain));
end
if( length(tree_node.left_constrain) > 0)
  y = y .* ((XData(tree_node.dim, :) > tree_node.left_constrain));
end
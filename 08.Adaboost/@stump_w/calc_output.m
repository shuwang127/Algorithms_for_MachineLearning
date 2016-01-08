%   This file is part of GML Matlab Toolbox
%   For conditions of distribution and use, see the accompanying License.txt file.

function y = calc_output(stump, XData)

y = (XData(stump.t_dim, :) <= stump.threshold) * (stump.signum) + (XData(stump.t_dim, :) > stump.threshold) * (-stump.signum);
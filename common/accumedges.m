function val4e = accumedges(e4ed, val4ed)
%%ACCUMEDGES accumulates global contributions on each edge for each 
%triangle. The data structures e4ed specifies the connectivity of edges 
%and triangles in a triangulation. The val4ed must be a vector with 
%the number of triangles in the last dimension. Values from interior edges
%are weighted by a factor 1/2.
%   val4e = ACCUMEDGES(e4ed, val4ed) creates array val4e (float: nElem x 1)
%   for the values in val4ed (float: nEdges x 1)
%

% Copyright 2021 Philipp Bringmann
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

    % Abbreviate geometric information
    nElem = max(e4ed(:));
    isInteriorEdge = (e4ed(:,2) ~= 0);

    % Accumulate values on interior edges
    val4e = accumarray(reshape(e4ed(isInteriorEdge,:), [], 1), ...
                       repmat(val4ed(isInteriorEdge) / 2, 2, 1), [nElem 1]);

    % Accumulate values on boundary edges
    val4e = val4e + accumarray(reshape(e4ed(~isInteriorEdge,1), [], 1), ...
                               val4ed(~isInteriorEdge), [nElem 1]);
end
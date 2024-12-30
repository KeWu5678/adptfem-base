function jump4ed = accumjumps(e4ed, ed4e, evaluation4e)
%%ACCUMJUMPS accumulates local contributions on each edge of a triangle
%into a global jump term. The data structures e4ed and ed4e specify the
%jumps in a triangulation. The evaluation4e must be an array with the 
%three dimensions: 1. number of point evaluations (e.g. quadrature points)
%on each edge; 2. one of the three edges; 3. number of triangles.
%   jump4ed = ACCUMJUMPS(e4ed, ed4e, evaluation4e) creates the
%   array jump4ed (float: nPoints x nEdges) for the values
%   in evaluation4e (float: nPoints x 3 x nElem)
%
%Warning 1:
%
%The contributions from both adjacent edges are ADDED. Make sure
%that the signs are appropriate to compute differences of values along the
%edge.
%
%Warning 2:
%
%Values of the point evaluations are added up in opposite order
%to match the same points in the counter-clockwise enumeration
%leading to opposite enumeration of points on both sides of the edge.
%This only works for symmetric points on the edges.
%
%Example:
%
%Given values v1, v2, v3 on edge E in T+ and values w1, w2, w3 on E in T-,
%the function returns the jump [v1 + w3, v2 + w2, v3 + w1] for the edge E
%as depicted below.
%
%        /\
%       /  \
%      /    \
%     /  T+  \
%    /        \
%   /          \
%  / v1  v2  v3 \
% o-----[E]------o
%  \ w3  w2  w1 /
%   \          /
%    \        /
%     \  T-  /
%      \    /
%       \  /
%        \/
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

    % Abbreviate evaluation information
    nPoints = size(evaluation4e, 1);

    % Abbreviate edge information
    nEdges = size(e4ed, 1);
    edgeNumbers = (1:nEdges)';
    isInterior = (e4ed(:,2) ~= 0);
    nInteriorEdges = nnz(isInterior);

    % Translate local into global edge numbers
    ed4ePlus = ed4e(e4ed(:,1),:);
    IndEdgePlus = permute(bsxfun(@eq, ed4ePlus, edgeNumbers), [3 2 1]);
    ed4eMinus = ed4e(e4ed(isInterior,2),:);
    IndEdgeMinus = ...
        permute(bsxfun(@eq, ed4eMinus, edgeNumbers(isInterior)), [3 2 1]);

    % Prepare translation of local into global edge numbers for specific
    % number of points in the evaluation
    IndEdgePlus = repmat(IndEdgePlus, [nPoints 1 1]);
    IndEdgeMinus = repmat(IndEdgeMinus, [nPoints 1 1]);

    % Extract evaluation on T_+
    evaluation4ePlus = evaluation4e(:,:,e4ed(:,1));
    evaluation4ePlus = ...
        reshape(evaluation4ePlus(IndEdgePlus), [nPoints nEdges]);

    % Extract evaluation on T_-
    evaluation4eMinus = evaluation4e(:,:,e4ed(isInterior,2));
    evaluation4eMinus = ...
        reshape(evaluation4eMinus(IndEdgeMinus), [nPoints nInteriorEdges]);

    % Accumulate jumps
    jump4ed = evaluation4ePlus;
    jump4ed(:,isInterior) = ...
        jump4ed(:,isInterior) + evaluation4eMinus(end:-1:1,:);
    
end
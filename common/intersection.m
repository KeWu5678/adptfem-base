function [A, polygon] = intersection(polygon1, polygon2)
%%INTERSECTION computes the area of intersection of the two convex polygons
%given by lists of vertices
%   [A, polygon] = INTERSECTION(polygon1, polygon2) returns the area A of
%   the intersection and its set of vertices

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

    %% DETERMINE VERTICES OF INTERSECTION POLYGON
    % by Sutherland–Hodgman algorithm

    % Initialisation
    nEdges = size(polygon1, 1);
    selected = polygon2;

    % Iterate over the edges of polygon1
    for j = 1:nEdges
        % The results from the previous iteration are the candidates for the
        % next edge
        candidates = selected;
        nCandidates = size(candidates, 1);
        selected = zeros(0, 2);

        % Current edge of polygon1
        edge1 = polygon1([j mod(j, nEdges)+1],:);

        % Iterate over vertices to check
        for k = 1:nCandidates
            % Extract current vertex and determine previous vertex in list of
            % candidates
            current = candidates(k,:);
            k_previous = mod(k-1, nCandidates);
            if k_previous == 0
                k_previous = nCandidates;
            end
            previous = candidates(k_previous,:);

            % Edge formed by the current and previous vertex in the intersection
            edge2 = [current; previous];

            % Compute intersection of edge1 and edge2
            tmp1 = [edge1(1) * edge1(4) - edge1(2) * edge1(3); ...
                    edge2(1) * edge2(4) - edge2(2) * edge2(3)];
            tmp2 = [edge1(1) - edge1(2); edge2(1) - edge2(2)];
            tmp3 = [edge1(3) - edge1(4); edge2(3) - edge2(4)];
            intersection_node = [tmp1(1) * tmp2(2) - tmp1(2) * tmp2(1),....
                                 tmp1(1) * tmp3(2) - tmp1(2) * tmp3(1)] ...
                                 / (tmp2(1) * tmp3(2) - tmp2(2) * tmp3(1));

            % Check if vertices lie inside wrt to current hyperplane
            det1 = current(1) * (edge1(3) - edge1(4)) ...
                   + current(2) * (edge1(2) - edge1(1)) ...
                   + edge1(1) * edge1(4) - edge1(3) * edge1(2);
            current_inside = (det1 >= 0);
            det2 = previous(1) * (edge1(3) - edge1(4)) ...
                   + previous(2) * (edge1(2) - edge1(1)) ...
                   + edge1(1) * edge1(4) - edge1(3) * edge1(2);
            previous_inside = (det2 >= 0);

            % Decide whether current and/or intersection node belong to the
            % intersection polygon
            if current_inside
                if ~previous_inside
                    selected(end+1,:) = intersection_node; %#ok<*AGROW>
                end
                selected(end+1,:) = current;
            elseif previous_inside
                selected(end+1,:) = intersection_node;
            end
        end
    end

    % Process result
    polygon = selected;
    nSelected = size(selected, 1);

    %% COMPUTE AREA OF INTERSECTION POLYGON
    if nSelected >= 3
        % Evaluate formula for area of a polygon
        ind_shift = [2:nSelected, 1];
        A = abs(sum(selected(:,1) .* selected(ind_shift,2) ...
                    - selected(ind_shift,1) .* selected(:,2))) / 2;
    else
        A = 0;
    end

end

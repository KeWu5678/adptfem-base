function ax = plotCR(meshdata, u)
%%PLOTCR plots a piecewise P1 that is continuous at the edge midpoints
%(Crouzeix-Raviart finite element function) given by the coefficient
%vector u of the values in the midpoints with respect to the
%triangulation meshdata
%   ax = PLOTCR(meshdata, u)

% Copyright 2020 Philipp Bringmann
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

    % Create axes object
    ax = gca;

    % Clear current axes object
    cla(ax);

    if size(u, 2) == 1
        % Rearrange coordinates of triangle vertices
        c4e = meshdata.c4e;
        X = squeeze(c4e(1,:,:));
        Y = squeeze(c4e(2,:,:));
        % Translate values in the edge midpoints into values in the nodes
        u4e = u(meshdata.ed4e)';
        Z = (ones(3) - 2*eye(3)) * u4e;
        % Determine the colour by the value in the midpoint of the triangle
        C = sum(u4e, 1) / 3;
        % Omit black edges for large number of triangles
        if( meshdata.nElem > 2000 )
            patch(X, Y, Z, C, 'EdgeColor', 'none');
        else
            patch(X, Y, Z, C);
        end
        % Set up 3D graph plot
        view(-37.5, 30);
    elseif size(u, 2) == 2
        % Plot of vector field
        set(ax, 'XLim', 1.2*[min(meshdata.c4n(:,1)), max(meshdata.c4n(:,1))],...
                'YLim', 1.2*[min(meshdata.c4n(:,2)), max(meshdata.c4n(:,2))]);
        quiver2(meshdata.mid4ed(:,1), meshdata.mid4ed(:,2), ...
                u(:,1), u(:,2), 'n=', 0.1, 'w=', [1 1]);
    else
        error('Invalid size of coefficient vector');
    end
    
    % Draw title
    title(ax, {'CR function plot'; [num2str(meshdata.nEdges), ' edges']});
end

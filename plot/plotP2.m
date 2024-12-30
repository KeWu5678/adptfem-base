function ax = plotP2(meshdata, x4e)
%%PLOTP2 plots a piecewise P2 function given by the coefficient vector x4e
%on the triangulation meshdata
%   ax = PLOTP2(meshdata, x4e)

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

    % Create axes object
    ax = gca;

    % Clear current axes object
    cla(ax);

    % No of triangles in given triangulation
    nElemOld = meshdata.nElem;

    % Uniform red-refinement of triangulation
    meshdata = meshdata.refineUniformRed();

    % P1 plot of nodal interpolation on red-refinement
    if size(x4e, 3) == 1
        % Rearrange coordinates of triangle vertices
        c4e = meshdata.c4e;
        X = squeeze(c4e(1,:,:));
        Y = squeeze(c4e(2,:,:));
        % Translate values in the edge midpoints into values in the nodes
        Z = zeros(3, nElemOld, 4);
        Z(:,:,1) = x4e(:,[1 6 5])';
        Z(:,:,2) = x4e(:,[6 2 4])';
        Z(:,:,3) = x4e(:,[5 4 3])';
        Z(:,:,4) = x4e(:,[4 5 6])';
        Z = reshape(Z, 3, []);
        % Determine the colour by the value in the midpoint of the triangle
        C = sum(Z, 1) / 3;
        % Omit black edges for large number of triangles
        if( meshdata.nElem > 2000 )
            patch(ax, X, Y, Z, C, 'EdgeColor', 'none');
        else
            patch(ax, X, Y, Z, C);
        end
    else
        error('Invalid size of coefficient vector');
    end

    % Set up 3D graph plot
    view(ax, -37.5, 30);

    % Draw title
    title(ax, {'P2 function plot'; [num2str(meshdata.nNodes), ' nodes']});
end

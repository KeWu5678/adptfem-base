function ax = plotMorley(meshdata, u)
%%PLOTMORLEY plots a scalar piecewise quadratic Morley finite element\
%function given by the coefficient vector u on the triangulation meshdata
%   ax = PLOTMORLEY(meshdata, u)

% Copyright 2022 Philipp Bringmann
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

    % Uniform red-refinement of triangulation
%     meshdata = meshdata.refineUniformRed();

    % Rearrange coordinates of triangle vertices
    c4e = meshdata.c4e;
    X = squeeze(c4e(1,:,:));
    Y = squeeze(c4e(2,:,:));

    % Extract local coefficients on each triangle
    u4e = u(meshdata.n4e)';
    Z = u4e;

    % Determine the colour by the value in the midpoint of the triangle
    C = sum(u4e, 1) / 3;

    % S1 plot of triangle-wise nodal interpolation on red-refinement
    if( meshdata.nElem > 2000 )
        patch(X, Y, Z, C, 'EdgeColor', 'none');
    else
        % Omit black edges for large number of triangles
        patch(X, Y, Z, C);
    end

    % Set up 3D graph plot
    view(-37.5, 30);

    % Draw title
    title(ax, {'Morley function plot'; [num2str(meshdata.nNodes), ' nodes']});
end

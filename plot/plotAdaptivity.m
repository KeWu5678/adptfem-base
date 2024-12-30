function ax = plotAdaptivity(meshdata)
%%PLOTP0 plots the fineness in terms of the mesh-size of a triangulation
%given by a MeshData object
%   ax = PLOTADAPTIVITY(meshdata)

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

    % Get the coordinates for each node of each triangle.
    X1 = meshdata.c4n(meshdata.n4e(:,1),1);
    Y1 = meshdata.c4n(meshdata.n4e(:,1),2);
    X2 = meshdata.c4n(meshdata.n4e(:,2),1);
    Y2 = meshdata.c4n(meshdata.n4e(:,2),2);
    X3 = meshdata.c4n(meshdata.n4e(:,3),1);
    Y3 = meshdata.c4n(meshdata.n4e(:,3),2);
    X = [X1'; X2'; X3'];
    Y = [Y1'; Y2'; Y3'];

    % Prepare function information
    Z = permute(repmat((sqrt(meshdata.area4e)), 1, 3), [2 1]);

    % Draw P0 function h4e as patch plot
    ax = gca;
    if meshdata.nElem > 2000
        patch(ax, X, Y, Z, Z, 'EdgeColor', 'none');
    else
        patch(ax, X, Y, Z, Z);
    end
    
    % Choose logarithmic scale for color
    set(ax, 'ColorScale', 'log');

    % Include colorbar
    colorbar(ax);

    % Set titles of plot
    title(ax, {'Adaptivity plot of mesh-size h_T'; ...
               [num2str(meshdata.nElem), ' elements']});
end


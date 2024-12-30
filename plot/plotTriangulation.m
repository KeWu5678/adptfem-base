function ax = plotTriangulation(meshdata)
%%PLOTTRIANGULATION plots the triangulation given by a MeshData object
%   ax = PLOTTRIANGULATION(meshdata)

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

    % Get the coordinates for each node of each triangle.
    X1 = meshdata.c4n(meshdata.n4e(:,1),1);
    Y1 = meshdata.c4n(meshdata.n4e(:,1),2);
    X2 = meshdata.c4n(meshdata.n4e(:,2),1);
    Y2 = meshdata.c4n(meshdata.n4e(:,2),2);
    X3 = meshdata.c4n(meshdata.n4e(:,3),1);
    Y3 = meshdata.c4n(meshdata.n4e(:,3),2);
    X = [X1'; X2'; X3'];
    Y = [Y1'; Y2'; Y3'];

    % Draw triangulation as patch plot
    patch(ax, X, Y, 'white', 'EdgeColor', 'blue');

    % Set titles of plot
    title(ax, {'Mesh plot'; [num2str(meshdata.nNodes), ' nodes']});
    
    % Format plot
    axis(ax, 'equal');
    axis(ax, 'tight');
end


function ax = plotRT0pw(meshdata, p4e)
%%PLOTRT0PW plots the piecewise Raviart-Thomas function p4e on a triangulation
%given by a MeshData object
%   ax = PLOTRT0PW(meshdata, p4e)

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

    % Get the coordinates for the midpoints of each triangle
    X = meshdata.mid4e(:,1);
    Y = meshdata.mid4e(:,2);

    % Prepare function information
    prefactors4e = ...
        bsxfun(@rdivide, meshdata.length4e, meshdata.area4e);
    c4e = meshdata.c4e;
    coefficients4e = ...
        [-sum(prefactors4e .* p4e .* permute(c4e(1,:,:), [3 2 1]), 2),...
         -sum(prefactors4e .* p4e .* permute(c4e(2,:,:), [3 2 1]), 2),...
          sum(prefactors4e .* p4e, 2)];
    U = coefficients4e(:,1) + coefficients4e(:,3) .* meshdata.mid4e(:,1);
    V = coefficients4e(:,2) + coefficients4e(:,3) .* meshdata.mid4e(:,2);

    % Draw triangulation as patch plot
    quiver2(X, Y, U, V, 'n=', 0.1, 'w=', [1 1]);
    xlim(ax, 1.1 * [min(meshdata.c4n(:,1)), max(meshdata.c4n(:,1))]);
    ylim(ax, 1.1 * [min(meshdata.c4n(:,2)), max(meshdata.c4n(:,2))]);

    % Set titles of plot
    title(ax, {'RT0pw plot'; [num2str(meshdata.nElem), ' elements']});
end


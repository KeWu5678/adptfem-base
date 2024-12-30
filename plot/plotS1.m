function ax = plotS1(meshdata, x)
%%PLOTS1 plots a piecewise P1 and globally continuous function (Courant
%finite element function) given by the coefficient vector x on the
%triangulation meshdata
%   ax = PLOTS1(meshdata, x)

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

    if size(x, 2) == 1
        % Plot of scalar function
        pat = trisurf(meshdata.n4e, meshdata.c4n(:,1),...
                      meshdata.c4n(:,2), x);
        if(meshdata.nElem > 2000)
            set(pat, 'EdgeColor', 'none');
        end
    elseif size(x, 2) == 2
        % Plot of vector field
        set(ax, 'XLim', 1.2*[min(meshdata.c4n(:,1)), max(meshdata.c4n(:,1))],...
                'YLim', 1.2*[min(meshdata.c4n(:,2)), max(meshdata.c4n(:,2))]);
        quiver2(meshdata.c4n(:,1), meshdata.c4n(:,2), ...
                x(:,1), x(:,2), 'n=', 0.1, 'w', [1 1]);
    else
        error('Invalid size of coefficient vector');
    end

    % Draw title
    title(ax, {'S1 function plot'; [num2str(meshdata.nNodes), ' nodes']});
end

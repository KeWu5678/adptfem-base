function ax = plotConvergence(ndof4lvl, error4lvl, legendEntry)
%%PLOTCONVERGENCE plots the error value with respect to the numbers of degrees
%of freedom in a double-logarithmic scale (convergence history plot)
%   ax = PLOTCONVERGENCE(ndof4lvl, error4lvl, legendEntry)

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

    % Display error for missing legendEntry
    if nargin < 3
        error('Please specify legend entry for graph');
    end

    % Generate a sequence of N^3 colours for the graphs
    N = 3;
    [red, green, blue] = meshgrid(linspace(0, 1, N));
    colourOrder = [red(:), green(:), blue(:)];

    % Save the current number of lines in the plot
    % for choosing a different colour
    fig = gcf;
    nLines = get(fig, 'UserData');
    if ~isempty(nLines)
        nLines = nLines + 1;
    else
        nLines = 1;
    end
    set(fig, 'UserData', nLines);

    % Create plot
    ax = gca;
    loglog(ax, ndof4lvl, error4lvl, '-s', ...
           'Color', colourOrder(nLines, :), 'DisplayName', legendEntry);

    % Add title
    title(ax, 'Convergence history plot');

    % Add axes labels
    xlabel(ax, 'ndof');
    ylabel(ax, 'error');

    % Update legend
    legend(ax, 'show');

    % Add new line into the current figure when calling plotConvergence again
    hold(ax, 'on');
end

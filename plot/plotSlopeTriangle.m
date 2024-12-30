function ax = plotSlopeTriangle(logx, logy, slope, swap)
%%PLOTSLOPETRIANGLE plots specified slope triangles in the current axis.
%This requires to specify the decadic logarithm of x and y coordinate of
%the base point of the triangle and the slope. All arguments allow for a
%list of values.
%
%   ax = PLOTSLOPETRIANGLE(logx, logy, slope) plots a lower triangle with
%       base point as the lower left vertex
%   ax = PLOTSLOPETRIANGLE(logx, logy, slope, swap) allows to plot the
%       (swapped) upper triangle with base point as the upper right vertex
%       if the boolean argument swap is true

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

    % Proceed input
    nTriangles = length(logx);

    % Create lower triangles if swapped is not specified
    if nargin < 4
        swap = false(nTriangles, 1);
    end

    % Create coordinates of slope triangles
    x = zeros(3, nTriangles);
    y = zeros(3, nTriangles);
    x(:,~swap) = [10.^logx(~swap), 10.^(logx(~swap) + 1), 10.^logx(~swap)]';
    y(:,~swap) = [10.^logy(~swap), 10.^logy(~swap), ...
                  10.^(logy(~swap) + slope(~swap))]';
    x(:,swap) = [10.^logx(swap), 10.^(logx(swap) - 1), 10.^logx(swap)]';
    y(:,swap) = [10.^logy(swap), 10.^logy(swap), ...
                 10.^(logy(swap) - slope(swap))]';

    % Create colour gradient
    c = zeros(3, nTriangles);
    c(:,~swap) = repmat([0.9; 0.9; 1], 1, nnz(~swap));
    c(:,swap) = repmat([1; 1; 0.9], 1, nnz(swap));
    caxis([0 1]);
    colormap(ax, "gray");

    % Plot slope triangles
    patch(x, y, c);

    % Create coordinates for text annotations
    xText1 = zeros(nTriangles, 1);
    xText2 = zeros(nTriangles, 1);
    xText1(~swap) = 0.9 * 10.^(logx(~swap) + 0.5);
    yText1(~swap) = 0.8 * 10.^logy(~swap);
%     xText2(~swap) = 0.8.^length(num2str(slope(~swap))) .* 10.^logx(~swap);
    xText2(~swap) = 1.1 * 10.^logx(~swap);
    yText2(~swap) = 10.^(logy(~swap) + slope(~swap)/2);
    xText1(swap) = 0.9 * 10.^(logx(swap) - 0.5);
    yText1(swap) = 1.2 * 10.^logy(swap);
    xText2(swap) = 1.1 * 10.^logx(swap);
    yText2(swap) = 10.^(logy(swap) - slope(swap)/2);

    % Plot text annotations
    text(xText1, yText1, "1");
    text(xText2, yText2, num2str(slope));

    % Sort data in axis
    data = get(ax, 'Children');
    newLines = [];
    rest = [];
    for j = 1:length(data)
        if contains(class(data(j)), 'Line')
            newLines(end+1) = data(j); %#ok<AGROW>
        else
            rest(end+1) = data(j); %#ok<AGROW>
        end
    end
    set(ax, 'Children', [rest, newLines]);

    % Remove entries for slope triangles from legend
    lgd = legend();
    lgdStrings = lgd.String;
    legend(lgdStrings(1:end-1));
end

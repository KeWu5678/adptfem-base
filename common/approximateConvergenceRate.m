function rate = approximateConvergenceRate(S, field)
%%APPROXIMATECONVERGENCERATE computes an approximated convergence rate using
%weighted least-squares linear fitting
%   rate = APPROXIMATECONVERGENCERATE(S) computes the rate of convergence of the
%   value 'eta' in the structure array S.level
%
%   rate = APPROXIMATECONVERGENCERATE(S, 'res') compute rate of convergence of
%   the specified field in the structure array S.level

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


    %% PROCEED INPUT
    if nargin < 2; field = 'eta'; end

    %% EXTRACT DATA
    nLvl = numel(S.level);
    ndof4lvl = zeros(nLvl, 1);
    value4lvl = zeros(nLvl, 1);
    for j = 1:nLvl
        ndof4lvl(j) = S.level(j).ndof;
        value4lvl(j) = S.level(j).(field);
    end

    %% COMPUTE WEIGHTING
    % weights = linspace(0, 1, nLvl);
    % weights = 2.^(0:nLvl-1);
    weights = 1 ./ (nLvl:-1:1);

    %% COMPUTE APPROXIMATED RATE OF CONVERGENCE
    X = [log(ndof4lvl), ones(nLvl, 1)];
    Y = log(value4lvl);
    linearModel = lscov(X, Y, weights);
    rate = - linearModel(1);

end


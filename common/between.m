function isBetween = between(x, y, z, tol)
%%BETWEEN checks whether x lies between y and z in every component
%   isBetween = BETWEEN(x, y, z) returns true if x > y - 1e4 * eps 
%   and x < z + 1e4 * eps
%
%   isBetween = BETWEEN(x, y, z, tol) returns true if x > y - tol 
%   and x < z + tol

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

    %% PROCEED INPUT
    if nargin < 4
        tol = 1e4 * eps;
    end

    %% CHECK WITH GIVEN TOLERANCE
    isBetween = (y - x) < tol & (x - z) < tol;
end

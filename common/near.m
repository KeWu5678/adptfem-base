function isNear = near(x, y, tol)
%%NEAR checks whether x and y are close in every component
%   isNeur = NEAR(x, y) returns true if |x - y| < 1e4 * eps
%
%   isNear = NEAR(x, y, tol) returns true if |x - y| < tol

% Copyright 2019 Philipp Bringmann
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
    if nargin < 3
        tol = 1e4 * eps;
    end

    %% COMPUTE ERROR
    err = abs(x - y);

    %% CHECK WITH GIVEN TOLERANCE
    isNear = err < tol;

end

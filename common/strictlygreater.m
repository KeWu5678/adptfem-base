function isGreater = strictlygreater(x, y, tol)
%%STRICTLYGREATER checks whether x is strictly greater than y in every
%component
%   isGreater = STRICTLYGREATER(x, y) returns true if x > y + 1e4 * eps
%
%   isGreater = STRICTLYGREATER(x, y, tol) returns true if x > y + tol

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

    %% COMPUTE DIFFERENCE
    dif = x - y;

    %% CHECK WITH GIVEN TOLERANCE
    isGreater = dif > tol;

end
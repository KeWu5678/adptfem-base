function Y = aitken(X, dim)
%%AITKEN computes Aitken delta^2 extrapolation for series acceleration
%along the specified dimension. This reduces the specified dimension
%   Y = AITKEN(X, dim)

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


    % Default dimension 1 (columnwise)
    if nargin < 2
        if isrow(X)
            dim = 2;
        else
            dim = 1;
        end
    end

    % Check for sufficient values along the dimension dim
    assert(size(X, dim) >= 3, 'Not enough values for series acceleration');

    % Permute X such that dim becomes the first dimension
    order = [dim, setdiff(1:ndims(X), dim)];
    X = permute(X, order);

    % Apply Aitken delta^2 formula
    Y = X(3:end,:) - (X(3:end,:) - X(2:end-1,:)).^2 ./ ...
        ((X(3:end,:) - X(2:end-1,:)) - (X(2:end-1,:) - X(1:end-2,:)));

    % Reverse permutation to obtain original order
    Y = ipermute(Y, order);
end

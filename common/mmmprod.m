function val = mmmprod(A, S, T)
%%MMMPROD computes the triple matrix product S' * A * T for a collection
%of matrices A, B, and C
%   val = MMMPROD(A, S, T) returns the product S' * A * T
%   val = MMMPROD(A, S) returns the product S' * A * S
%
%Matrix | Dimension
%------------------
%   A   | K x L x N
%   S   | K x J x N
%   T   | L x M x N
%  val  | J x M x N
%
%The third dimension may be expanded if it is a singleton. 

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

    % Proceed input
    if nargin < 3
        T = S;
    end

    % Check input data
    assert(size(A, 1) == size(S, 1), ...
           "Dimensions of A and S are not compatible");
    assert(size(A, 2) == size(T, 1), ...
           "Dimensions of A and T are not compatible");
    assert(size(A, 3) == 1 | size(S, 3) == 1 | size(A, 3) == size(S, 3), ...
           "Numbers of matrices in A and S are not compatible");
    assert(size(A, 3) == 1 | size(T, 3) == 1 | size(A, 3) == size(T, 3), ...
           "Numbers of matrices in A and T are not compatible");
    assert(size(S, 3) == 1 | size(T, 3) == 1 | size(S, 3) == size(T, 3), ...
           "Numbers of matrices in S and T are not compatible");

    % Computation
    val = sum(bsxfun(@times, permute(A, [1 2 4 3]), permute(T, [4 1 2 3])), 2);
    val = sum(bsxfun(@times, permute(S, [1 2 4 3]), val), 1);
    val = squeeze(val);

end

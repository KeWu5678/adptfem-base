function val = vmvprod(A, u, v)
%%VMVPROD computes the vector-matrix-vector product u' * A * v
%for a collection of M vectors u,v and a collection of M x M matrix A
%   val = vmvprod(A, u, v) returns an N vector for given M x M x N array A
%   and N x M arrays u and v

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

    val = bsxfun(@times, permute(A, [3 1 2]), permute(v, [1 3 2]));
    val = squeeze(sum(val, 3));
    val = sum(val.*u, 2);

end

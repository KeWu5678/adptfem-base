function A = assemble(A4p, dof4p, N)
%%ASSEMBLE carries out the assembling of global matrix A from local
%contributions from A4p. The relation between local and global indices is
%described by the array dof4p. In order to save memory, the assembling is 
%carried out in serial and all GPU arrays are gathered first.
%
%   A = ASSEMBLE(A4p, dof4p, N) creates a quadratic N-by-N matrix A
%   A = ASSEMBLE(A4p, dof4p) determines N as the maximal index in dof4p

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

    % Gather possibly distributed arrays
    [A4p, dof4p] = gather(A4p, dof4p);

    % Compute number of local degrees of freedom
    ndof4p = size(dof4p, 1);

    % Optionally compute number dimension of output array
    if nargin < 3
        N = max(dof4p(:));
    end

    % Assembling
    if ismatrix(A4p)
        A = accumarray(dof4p(:), A4p(:), [N 1]);
    elseif ndims(A4p) == 3
        RowIndices4p = permute(repmat(dof4p, 1, 1, ndof4p), [1 3 2]);
        ColumnIndices4p = permute(RowIndices4p, [2 1 3]);
        A = sparse(RowIndices4p(:), ColumnIndices4p(:), A4p(:), N, N);
    else
        error('Multidimensional assembling not supported');
    end
end
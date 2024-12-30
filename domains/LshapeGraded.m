function [c4n, n4e, onDirichlet, onNeumann] = LshapeGraded(beta, N)
%%LSHAPEGRADED returns the triangulation of the L-shaped domain with
%beta-grading
%   [c4n, n4e, onDirichlet, onNeumann] = LshapeGraded(beta, N)

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

    % Define macro triangulation
    c4n_ma = [-1 -1; 0 -1; -1 0; 0 0; 1 0; -1 1; 0 1; 1 1];
    n4e_ma = [4 5 8; 4 8 7; 4 7 6; 4 6 3; 4 3 1; 4 1 2];

    % Transform graded triangulation to the six triangles in the L-shaped
    % domain
    [c4n_mi, n4e_mi, nE1, ~, nE3] = geomTrefG(N, beta);
    L = size(c4n_mi,1); n4e = []; c4n = [];
    for j = 1:6
        c4n_tmp = (1-c4n_mi(:,1)-c4n_mi(:,2)) * c4n_ma(n4e_ma(j,1),:) ...
                  + c4n_mi(:,1) * c4n_ma(n4e_ma(j,2),:) ...
                  + c4n_mi(:,2) * c4n_ma(n4e_ma(j,3),:);
        n4e = [n4e; n4e_mi+(j-1)*L];
        c4n = [c4n; c4n_tmp];
    end

    % Remove duplettes
    idn = [ones(1,5); (L+1):L:(5*L+1)]';
    nE1(1,:) = []; nE3(1,:) = [];
    for j = 1:5
        idn = [idn;[unique(nE3(:))+(j-1)*L, unique(nE1(:))+j*L]];
    end
    for j = 1:size(idn, 1)
        n4e(n4e == idn(j,2)) = idn(j,1);
    end
    ctr = 0;
    for j = 1:size(c4n,1)
        ind = (n4e == j);
        if any(ind(:))
            ctr = ctr + 1;
            n4e(ind) = ctr;
        end
    end
    c4n(idn(:,2),:) = [];

    % Transform domain from Gamma- to L-shaped domain
    c4n = [-c4n(:,2), c4n(:,1)];

    % Boundary information
    onDirichlet =  @onDirichletBoundary;

    onNeumann = @(x) false(size(x, 1), 1);
end


function [c4n, n4e, nE1, nE2, nE3] = geomTrefG(N, beta)
%%GEOMTREFG creates graded triangulation of reference triangle
    c4n = [0,0]; n4e = [1 2 3];
    for j = 1 : N
        xi = (j/N)^beta;
        for m = 0 : j
            c4n = [c4n;[xi,0]+(m/j)*[-xi,xi]];
        end
    end
    for j = 1 : N-1
        for n = 1 : j
            n4e = [n4e;j*(j+1)/2+n+[0,j+2,1]];
        end
        for n = 1 : j+1
            n4e = [n4e;j*(j+1)/2+n+[0,j+1,j+2]];
        end
    end
    nE1 = [(0:N-1).*(1:N)/2+1;(1:N).*(2:N+1)/2+1]';
    nE2 = N*(N+1)/2+[1:N; 2:N+1]';
    nE3 = [(1:N).*(2:N+1)/2;(2:N+1).*(3:N+2)/2]';
end


function onDb = onDirichletBoundary(x)
%%ONDIRICHLETBOUNDARY checks whether given points belong to the Dirichlet
%boundary
    onDb = near(x(:,1), -1) | near(x(:,1), 1) | ...
           near(x(:,2), -1) | near(x(:,2), 1) | ...
           (near(x(:,1), 0) & nearorgreater(x(:,2), 0)) | ...
           (near(x(:,2), 0) & nearorgreater(x(:,1), 0));
end
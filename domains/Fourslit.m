function [c4n, n4e, onDirichlet, onNeumann] = Fourslit(reflected)
%%FOURSLIT provides geometric information of four slit domain.
%According to the optional boolean argument reflected, the triangulation may
%be reflected or follow legacy counter-clockwise numbering conventions with
%refinement edge between first two local nodes.
%   [c4n, n4e, onDirichlet, onNeumann] = FOURSLIT(reflected)

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
% along with this prograGa  If not, see <http://www.gnu.org/licenses/>.
%


    % Proceed input
    if nargin == 0
        reflected = false;
    end

    % Coordinates of nodes
    [Y, X] = meshgrid(linspace(-1, 1, 5), linspace(-1, 1, 5));
    c4n = [[X(:), Y(:)];
            0, -1;
           -1,  0;
            1,  0;
            0,  1];

    % Node numbers of the elements
    if reflected
        % reflected triangulation
        n4e = [ 6  2  1;
                6  2  7;
                8  2  7;
                8  2  3;
                8  4 26;
                8  4  9;
               10  4  9;
               10  4  5;
                6 12 11;
                6 12  7;
                8 12  7;
                8 12 13;
                8 14 13;
                8 14  9;
               10 14  9;
               10 14 15;
               16 12 27;
               16 12 17;
               18 12 17;
               18 12 13;
               18 14 13;
               18 14 19;
               20 14 19;
               20 14 28;
               16 22 21;
               16 22 17;
               18 22 17;
               18 22 23;
               18 24 29;
               18 24 19;
               20 24 19;
               20 24 25];
    else
        % triangles in counter-clockwise enumeration
        n4e = [ 2  6  1;
                6  2  7;
                2  8  7;
                8  2  3;
                4  8 26;
                8  4  9;
                4 10  9;
               10  4  5;
                6 12 11;
               12  6  7;
                8 12  7;
               12  8 13;
                8 14 13;
               14  8  9;
               10 14  9;
               14 10 15;
               12 16 27;
               16 12 17;
               12 18 17;
               18 12 13;
               14 18 13;
               18 14 19;
               14 20 19;
               20 14 28;
               16 22 21;
               22 16 17;
               18 22 17;
               22 18 23;
               18 24 29;
               24 18 19;
               20 24 19;
               24 20 25];
    end

    % Boundary information
    onDirichlet =  @onDirichletBoundary;

    onNeumann = @(x) false(size(x, 1), 1);

end


function onDb = onDirichletBoundary(x)
%%ONDIRICHLETBOUNDARY checks whether given points belong to the Dirichlet
%boundary
    onDb = near(x(:,1), -1) | near(x(:,1), 1) | ...
           near(x(:,2), -1) | near(x(:,2), 1) | ...
           (near(x(:,1), 0) & ...
              (nearorless(x(:,2), -0.5) | nearorgreater(x(:,2), 0.5))) | ...
           (near(x(:,2), 0) & ...
              (nearorless(x(:,1), -0.5) | nearorgreater(x(:,1), 0.5)));
end

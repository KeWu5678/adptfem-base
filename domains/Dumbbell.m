function [c4n, n4e, onDirichlet, onNeumann] = Dumbbell(reflected)
%%DUMBBELL provides geometric information of 1/8 cusp domain.
%According to the optional boolean argument reflected, the triangulation may
%be reflected or follow legacy counter-clockwise numbering conventions with
%refinement edge between first two local nodes.
%   [c4n, n4e, onDirichlet, onNeumann] = DUMBBELL(reflected)

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
    c4n = [ -1     -1;
            -1      0;
            -1      0;
            -1      1;
             0     -1;
             0      0;
             0      1;
             0.5   -1;
             0.5   -0.5;
             0.75  -1;
             0.75  -0.75;
             1     -1;
             1     -0.75;
             1     -0.5;
             1      0;
             1      1;
             1.25  -1;
             1.25  -0.75;
             1.5   -1;
             1.5   -0.75;
             1.75  -1;
             1.75  -0.75;
             2     -1;
             2     -0.75;
             2.25  -1;
             2.25  -0.75;
             2.5   -1;
             2.5   -0.75;
             2.75  -1;
             2.75  -0.75;
             3     -1;
             3     -0.75;
             3     -0.5;
             3      0;
             3      1;
             3.25  -1;
             3.25  -0.75;
             3.5   -1;
             3.5   -0.5;
             4     -1;
             4      0;
             4      1;
             5     -1;
             5      0;
             5      1];

    % Node numbers of the elements
    if reflected
        % reflected triangulation
        n4e = [ 1  5  6;
                1  2  6;
                4  3  6;
                4  7  6;
               16  7  6;
               16 15  6;
                9 15  6;
                9  5  6;
                9  5  8;
                9 15 14;
                9 11  8;
                9 11 14;
               10 11  8;
               13 11 14;
               10 11 12;
               13 11 12;
               13 17 12;
               13 17 18;
               19 17 18;
               19 20 18;
               19 20 21;
               22 20 21;
               22 23 21;
               22 23 24;
               26 23 24;
               26 23 25;
               26 28 25;
               27 28 25;
               27 28 30;
               27 29 30;
               32 29 30;
               32 29 31;
               32 37 31;
               36 37 31;
               32 37 33;
               36 37 38;
               39 37 33;
               39 37 38;
               39 34 33;
               39 40 38;
               39 30 41;
               39 34 41;
               35 34 41;
               35 42 41;
               45 42 41;
               45 44 41;
               43 44 41;
               43 40 41];
    else
        % triangles in counter-clockwise enumeration
        n4e = [ 6  1  5;
                1  6  2;
                6  4  3;
                4  6  7;
                6 16  7;
               16  6 15;
               15  6  9;
                6  5  9;
                9  5  8;
               15  9 14;
                9  8 11;
               14  9 11;
               11  8 10;
               14 11 13;
               12 11 10;
               11 12 13;
               17 13 12;
               13 17 18;
               19 18 17;
               18 19 20;
               21 20 19;
               20 21 22;
               23 22 21;
               22 23 24;
               23 26 24;
               26 23 25;
               25 28 26;
               28 25 27;
               27 30 28;
               30 27 29;
               29 32 30;
               32 29 31;
               31 37 32;
               37 31 36;
               37 33 32;
               38 37 36;
               39 33 37;
               38 39 37;
               39 34 33;
               40 39 38;
               40 41 39;
               41 34 39;
               41 35 34;
               35 41 42;
               41 45 42;
               45 41 44;
               41 43 44;
               43 41 40];
    end

    % Boundary information
    onDirichlet =  @onDirichletBoundary;

    onNeumann = @(x) false(size(x, 1), 1);

end


function onDb = onDirichletBoundary(x)
%%ONDIRICHLETBOUNDARY checks whether given points belong to the Dirichlet
%boundary
    onDb = near(x(:,1), -1) | near(x(:,1), 5) | ...
           near(x(:,2), -1) | near(x(:,2), 1) | ...
           (near(x(:,1), 1) & nearorgreater(x(:,2), -0.75)) | ...
           (near(x(:,1), 3) & nearorgreater(x(:,2), -0.75)) | ...
           (nearorgreater(x(:,1), 1) & nearorless(x(:,1), 3) ...
            & near(x(:,2), -0.75));
end


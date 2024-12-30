function [c4n, n4e, onDirichlet, onNeumann, normal] = Rectangle(reflected)
%%RECTANGLE provides the rectangle geometry with pure Dirichlet boundary.
%According to the optional boolean argument reflected, the triangulation may
%be reflected or follow legacy counter-clockwise numbering conventions with
%refinement edge between first two local nodes.
%   [c4n, n4e, onDirichlet, onNeumann, normal] = RECTANGLE(reflected)

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
% along with this prograGa  If not, see <http://www.gnu.org/licenses/>.
%


    % Proceed input0
    if nargin == 0
        reflected = false;
    end

    % DOMAIN
    % Coordinates of nodes
    c4n = [-1   0;
            0   0;
            1   0;
           -0.5 0.5;
            0.5 0.5;
           -1   1;
            0   1;
            1   1];

    % Node numbers of the elements
    if reflected
        % reflected triangulation
        n4e = [ 1   4   2;
                1   4   6;
                7   4   6;
                7   4   2;
                7   5   2;
                7   5   8;
                3   5   8;
                3   5   2];
    else
        % triangles in counter-clockwise enumeration
        n4e = [ 1   2   4;
                6   1   4;
                7   6   4;
                2   7   4;
                7   2   5;
                8   7   5;
                3   8   5;
                2   3   5];
    end

    % Boundary information
    onDirichlet =  @onDirichletBoundary;

    onNeumann = @(x) false(size(x, 1), 1);

    normal = @normalVector;

end


function onDb = onDirichletBoundary(x)
%%ONDIRICHLETBOUNDARY checks whether given points belong to the Dirichlet
%boundary
    onDb = near(x(:,1), -1) | near(x(:,1), 1) | ...
           near(x(:,2), 0) | near(x(:,2), 1);
end


function normal = normalVector(x)
%%NORMALVECTOR computes the outward unit normal in x
%This function requires x to be an element of the
%closure of the domain Omega
    normal = zeros(size(x, 1), 2);
    normal(near(x(:,1), 0),1) = -1;
    normal(near(x(:,2), 0),2) = -1;
    normal(near(x(:,1), 1),1) = 1;
    normal(near(x(:,2), 1),2) = 1;
end

function [c4n, n4e, onDirichlet, onNeumann, normal] = Triangle(c4n)
%%TRIANGLE provides the triangle geometry with pure Dirichlet boundary.
%The optional argument c4n allows to create a triangulation of an arbitrary
%triangle instead of the default reference triangle.
%   [c4n, n4e, onDirichlet, onNeumann, normal] = TRIANGLE(reflected, c4n)

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
        c4n = [ 1    0;
                0    1;
                0    0];
    end

    % Node numbers of the elements
    n4e = [ 2   3   1];

    % Boundary information
    onDirichlet =  @(x) onDirichletBoundary(x, c4n);

    onNeumann = @(x) false(size(x, 1), 1);

    normal = @(x) normalVector(x, c4n);

end


function onDb = onDirichletBoundary(x, c4n)
%%ONDIRICHLETBOUNDARY checks whether given points belong to the Dirichlet
%boundary

    % Transform every point back to reference triangle
    transformation = [c4n(1,:) - c4n(3,:); c4n(2,:) - c4n(3,:)];
    xref = (x - c4n(3,:)) / transformation;
    % Check position with respect to reference triangle
    onDb = near(xref(:,1), 0) | near(xref(:,2), 0) | ...
           near(xref(:,1) + xref(:,2), 1);
end


function normal = normalVector(x, c4n)
%%NORMALVECTOR computes the outward unit normal in x
%This function requires x to be an element of the
%closure of the domain Omega

    % Compute normals of the three edges
    tangents = [c4n([2 3 1],:) - c4n];
    normals = [tangents(:,2), -tangents(:,1)];
    % Transform every point back to reference triangle
    transformation = [c4n(1,:) - c4n(3,:); c4n(2,:) - c4n(3,:)];
    xref = (x - c4n(3,:)) / transformation;
    % Assign corresponding normals
    normal = zeros(size(x, 1), 2);
    normal(near(xref(:,2), 0),:) = normals(1,:);
    normal(near(xref(:,1) + xref(:,2), 1),:) = normals(2,:);
    normal(near(xref(:,1), 0),:) = normals(3,:);
end

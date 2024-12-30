classdef Quadrature < handle
%%QUADRATURE objects provide quadrature points and weights for numerical
%integration on the reference domain in 1D (interval), 2D (triangle), or 3D
%(tetrahedron). The rules are exact at least for polynomials of specified
%partial degree in spatial dimension dim
%
%   q = QUADRATURE(dim, degree)
%   q.points
%   q.weights

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


    properties
        % spatial dimension (int: 1)
        dim
        % accuracy up to polynomial degree (int: 1)
        degree
        % quadrature points (float: nPoints x dim)
        points
        % quadrature weights (float: nPoints)
        weights
        % quadrature points transformed to each element (float: dim x nElem x nPoints)
        qp4e = nan
        % quadrature points transformed to each edge (float: dim x nEdges x nPoints)
        qp4ed = nan
    end

    methods
        %% CONSTRUCTOR
        function obj = Quadrature(dim, degree)
            %%QUADRATURE creates a Quadrature object
            %   q = QUADRATURE(dim, degree)
            obj.dim = dim;
            obj.degree = degree;
            obj.computePointsAndWeights();
        end

        %% BASIC FUNCTIONS
        function n = nPoints(obj)
            %%NPOINTS returns number of quadrature points
            n = size(obj.points, 1);
        end

        function clearMesh(obj)
            %%CLEARMESH deletes all quadrature nodes transformed
            %to a specific mesh
            obj.qp4e = nan;
            obj.qp4ed = nan;
        end

        %% INTEGRATION
        function int4p = integrate(obj, integrand4p, meas4p)
            %%INTEGRATE wrapper for quadrature evaluation
            %   int4p = QUADRATURE.INTEGRATE(integrand4p, meas4p) computes the
            %   integral for the already evaluated integrand4p
            weights4p = bsxfun(@times, permute(obj.weights, [2 3 4 5 1]),...
                              permute(meas4p, [2 3 4 1 5]));
            integrand4p = sum(sum(integrand4p, 2), 1);
            int4p = sum(weights4p .* integrand4p, 5);
            int4p = squeeze(int4p);
        end

        %% TRANSFORM QUADRATURE POINTS TO MESH
        function qp4e = getQp4e(obj, meshdata)
            %%GETQP4E transforms quadrature points from the reference
            %simplex to every simplex in the triangulation
            %   qp4e = QUADRATURE.GETQP4E(meshdata)
            %
            % Dimensions: 2 x nElem x nPoints
            if isnan(obj.qp4e)
                obj.qp4e = meshdata.transform(obj.points);
            end
            qp4e = obj.qp4e;
        end

        function qp4ed = getQp4ed(obj, meshdata)
            %%GETQP4ED transforms quadrature points from the reference
            %interval to every edge in the triangulation
            %   qp4ed = QUADRATURE.GETQP4ED(meshdata)
            %
            % Dimensions: 2 x nEdges x nPoints
            if obj.dim ~= 1
                error('Edge quadrature points not avaible by this object')
            end
            if isnan(obj.qp4ed)
                obj.qp4ed = ...
                    meshdata.edgeTransform(obj.points, 1:meshdata.nEdges);
            end
            qp4ed = obj.qp4ed;
        end

        function [qpPlus, qpMinus] = getQp4jumps(obj, meshdata)
            %%QP4JUMPS transforms quadrature points from the reference interval
            %to the reference triangle with respect to the adjacent simplices of
            %each face
            %   [qpPlus, qpMinus] = QUADRATURE.QP4JUMPS(meshdata)
            %
            % Dimensions: 2 x nFaces/nEdges x nPoints

            %% EDGE JUMPS
            if obj.dim == 1
                % Mesh information
                isInterior = (meshdata.e4ed(:,2) ~= 0);
                qp4ed = obj.getQp4ed(meshdata); %#ok<*PROPLC>

                % Transform points on edges back to reference triangle
                qpPlus = meshdata.inverseTransform(qp4ed, meshdata.e4ed(:,1));
                qpMinus = zeros(size(qpPlus));
                qpMinus(:,isInterior,:) = ...
                    meshdata.inverseTransform(qp4ed(:,isInterior,:), ...
                                              meshdata.e4ed(isInterior,2));
            else
                error(['Quadrature.qp4jumps not implemented for dimension',...
                       num2str(obj.dim)]);
            end
        end

        %% EVALUATIONS OF FUNCTION HANDLES IN GAUSS POINTS
        function val4e = evaluateFunction(obj, meshdata, fun)
            %%EVALUATEFUNCTION evaluates a function handle in every quadrature
            %point on every simplex
            %   val4e = QUADRATURE.EVALUATEFUNCTION(meshdata, fun)

            % Evaluate function
            qp4e = obj.getQp4e(meshdata);
            qp4e = permute(qp4e(:,:), [2 1]);
            val4e = fun(qp4e);
            % Reshape according to shape of the function
            comps = size(fun(qp4e(1,:)));
            if comps(1) == 1
                val4e = permute(val4e, [2 1]);
                comps = fliplr(comps);
            end
            val4e = reshape(val4e, [comps 1 meshdata.nElem obj.nPoints]);
        end

        function val4ed = edgeEvaluateFunction(obj, meshdata, fun, edges)
            %%EDGEEVALUATEFUNCTION evaluates a function handle in every
            %quadrature point on the given edges
            %   val4ed = QUADRATURE.EDGEEVALUATEFUNCTION(meshdata, fun, edges)

            % Evaluate function
            qp4ed = obj.getQp4ed(meshdata);
            qp4ed = qp4ed(:,edges,:);
            qp4ed = permute(qp4ed(:,:), [2 1]);
            val4ed = fun(qp4ed);

            % Reshape according to shape of the function
            comps = size(fun(qp4ed(1,:)));
            if comps(1) == 1
                val4ed = permute(val4ed, [2 1]);
                comps = fliplr(comps);
            end
            val4ed = reshape(val4ed, [comps 1 length(edges) obj.nPoints]);
        end
    end

    methods (Access = private)
        function computePointsAndWeights(obj)
            %%COMPUTEPOINTSANDWEIGHTS
            nPoints = ceil((obj.degree + 1) / 2);
            nPoints = max(nPoints, 1);  % compute at least one quadrature point

            %% COMPUTE GAUSS POINTS
            if obj.dim == 1
                [obj.points, obj.weights] = obj.computeGaussJacobi(nPoints, 0);
            elseif obj.dim == 2
                % Gauss-Jacobi points for unit interval
                [x1, w1] = obj.computeGaussJacobi(nPoints, 1);
                [x2, w2] = obj.computeGaussJacobi(nPoints, 0);

                % Tensor product Gauss points for unit square
                x = reshape(repmat(x1', nPoints, 1), [], 1);
                y = repmat(x2, nPoints, 1);

                % Gauss points and weights for reference triangle
                obj.points = [x, (1 - x) .* y];
                obj.weights = kron(w1, w2);

                % Divide by volume of reference triangle
                obj.weights = 2 * obj.weights;
            elseif obj.dim == 3
                % Gauss-Jacobi points for unit interval
                [x1, w1] = obj.computeGaussJacobi(nPoints, 2);
                [x2, w2] = obj.computeGaussJacobi(nPoints, 1);
                [x3, w3] = obj.computeGaussJacobi(nPoints, 0);

                % Gauss points for unit cube
                x = reshape(repmat(x1', nPoints^2, 1), [], 1);
                y = reshape(repmat(x2', nPoints, nPoints), [], 1);
                z = repmat(x3, nPoints^2, 1);

                % Gauss points for reference tetrahedron
                obj.points = [x, (1 - x) .* y, (1 - x) .* (1 - y) .* z];
                obj.weights = kron(w1, kron(w2, w3));

                % Divide by volume of reference triangle
                obj.weights = 6 * obj.weights;
            else
                error('Quadrature: Invalid dimension');
            end
        end
    end

    methods (Static, Access = private)
        function [x, w] = computeGaussJacobi(n, a)
            %%COMPUTEGAUSSJACOBI gives n Gauss-Jacobi points for the unit
            %interval with respect to the weight function w(x) = (1 - x)^a
            %computed by Golub-Welsh method
            %
            %   [x, w] = COMPUTEGAUSSJACOBI(n, a)

            % Recurrence terms in the formula
            %   p_j = (A_j x + B_j) p_j-1 - C_j p_j-2
            j = (1:n);
            D = 2 * j .* (j + a) .* (2*j + a - 2);
            D(D == 0) = 1;
            A = (2 * j + a - 2) .* (2 * j + a - 1) .* (2 * j + a) ./ D;
            A(A == 0) = 1;
            B = a^2 * (2 * j + a - 1) ./ D;
            C = 2 * (j + a - 1) .* (j - 1) .* (2 * j + a) ./ D;

            % Compute recurrence coefficients gamma and delta
            % by formulas from [Golub-Welsh]
            delta = - B ./ A;
            gamma = sqrt(C(2:n) ./ A(1:n-1) ./ A(2:n));

            % Find n Gauss-Jacobi points for interval [-1,1]
            J = diag(gamma, -1) + diag(delta) + diag(gamma, 1);
            [V, D] = eig(J);
            x = diag(D);
            w = 2^(a + 1) / (a + 1) * V(1,:).^2;

            % Linear map onto interval [0,1] and modification of weights
            x = 0.5 * x + 0.5;
            w = 0.5^(a + 1) * w';
        end
    end
end

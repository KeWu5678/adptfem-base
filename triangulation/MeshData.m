classdef MeshData < handle
%%MESHDATA class reprenting the geometric information of a mesh

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


    properties
        % coordinates of the nodes (float: nNodes x 2)
        c4n
        % indices of the 3 nodes of each triangle (int: nElem x 3)
        n4e
        % function handle to determine the location of coordinates
        onDirichlet
        % function handle to determine the location of coordinates
        onNeumann
        % area of each triangle (float: nElem x 1)
        area4e = nan
        % interior angles of each triangle (float: nElem x 3)
        angle4e = nan
        % length of each edges (float: nEdges x 1)
        length4ed = nan
        % mesh size in terms of area of the adjacent triangle T+ (float: nEdges x 1)
        h4ed = nan
        % indices of (at most) 3 neighbouring triangles (int: nElem x 3)
        neighbours4e = nan
        % affine transformation matrix for each triangle (float: 2 x 2 x nElem)
        Trafo4e = nan
        % inverse of affine transformation matrix (float: 2 x 2 x nElem)
        InvTrafo4e = nan
        % affine transformation matrix for each edge (float: 2 x 1 x nEdges)
        Trafo4ed = nan
        % barycenter of each triangle (float: nElem x 2)
        mid4e = nan
        % midpoint of each edge (float: nEdges x 2)
        mid4ed = nan
        % the 3 outward unit normal vectors of each triangle (float: 2 x 3 x nElem)
        normal4e = nan
        % the 3 unit tangent vectors along the edges of each triangle (float: 3 x 2 x nElem)
        tangent4e = nan
        % unit normal vector of each edge (float: nEdges x 2)
        normal4ed = nan
        % unit tangential vector of each edge (float: nEdges x 2)
        tangent4ed = nan
        % sign change of normal vectors of edges and triangles (float: nElem x 3)
        sign4e = nan
        % sign change of normal vectors of Dirichlet edges and triangles (float: nDbEdges x 1)
        sign4edDb = nan
        % sign change of normal vectors of Neumann edges and triangles (float: nNbEdges x 1)
        sign4edNb = nan
        % indices of the 2 nodes of each edge (int: nEdges x 2)
        n4ed = nan
        % indices of the 3 edges of each triangle (int: nElem x 3)
        ed4e = nan
        % indices of adjacent triangles for each edge (int: nEdges x 2)
        e4ed = nan
        % indices of three nodes in T+ and T- (int: nEdges x 3 x 2)
        n4edgePatch = nan
        % indices of three edges in T+ and T- (int: nEdges x 3 x 2)
        ed4edgePatch = nan
        % indices of nodes opposite to each edge on T+ and T- (int: nEdges x 2)
        oppositeNode4ed = nan
        % indices of two edges opposite to each edge on T+ and T- (int: nEdges x 2 x 2)
        oppositeEdges4ed = nan
        % indices of the nodes on the Dirichlet boundary (int: nDbNodes x 1)
        DbNodes = nan
        % indices of the Dirichlet edges (int: nDbEdges x 1)
        DbEdges = nan
        % indices of the nodes on the Neumann boundary (int: nNbNodes x 1)
        NbNodes = nan
        % indices of the Neumann edges (int: nNbEdges x 1)
        NbEdges = nan
        % indices of parent triangle (int: nElem x 1)
        parents4e = nan
        % local children index in parent triangle (int: nElem x 1)
        childNo4e = nan
        % indices of child triangles in refinement (int: nElem x 4)
        child4e = nan
        % number of children existing in refinement (int: nElem x 1)
        nChildren4e = nan
    end

    properties (Access = private)
        % indices of two nodes of each edge sorted lexicographically (int: nEdges x 2)
        n4edSorted = nan
    end

    methods
        %% CONSTRUCTOR
        function obj = MeshData(c4n, n4e, onDirichlet, onNeumann, ...
                                parents4e, childNo4e)
            %%MESHDATA constructs a MeshData object
            %   M = MESHDATA(c4n, n4e, onDirichlet, onNeumann)
            %   M = MESHDATA(c4n, n4e, onDirichlet, onNeumann, ...
            %                parents4e, childNo4e)
            if nargin == 4
                obj.c4n = c4n;
                obj.n4e = n4e;
                obj.onDirichlet = onDirichlet;
                obj.onNeumann = onNeumann;
            elseif nargin == 6
                obj.c4n = c4n;
                obj.n4e = n4e;
                obj.onDirichlet = onDirichlet;
                obj.onNeumann = onNeumann;
                obj.parents4e = parents4e;
                obj.childNo4e = childNo4e;
            else
                error('Invalid MeshData constructor call');
            end
        end

        %% GET FUNCTIONS FOR NUMBERS OF GEOMETRIC ENTITIES
        function n = nElem(obj)
            %%NELEM returns total number of triangles in the triangulation
            n = size(obj.n4e, 1);
        end

        function n = nEdges(obj)
            %%NEDGES returns total number of edges in the triangulation
            n = size(obj.n4ed, 1);
        end

        function n = nDbEdges(obj)
            %%NDBEDGES returns total number of Dirichlet edges
            n = size(obj.n4edDb, 1);
        end

        function n = nNbEdges(obj)
            %%NNBEDGES returns total number of Neumann edges
            n = size(obj.n4edNb, 1);
        end

        function n = nNodes(obj)
            %%NNODES returns total number of nodes in the triangulation
            n = size(obj.c4n, 1);
        end

        %% GET FUNCTIONS FOR GEOMETRIC MEASURES
        function area4e = get.area4e(obj)
            if isnan(obj.area4e)
                % Compute the area by the determinant formula
                c4e = obj.c4e;
                obj.area4e = squeeze(abs(...
                    sum((c4e(2,[2 3 1],:) - c4e(2,[3 1 2],:)) ...
                        .* c4e(1,:,:), 2))) / 2;
            end
            area4e = obj.area4e;
        end

        function angle4e = get.angle4e(obj)
            if isnan(obj.angle4e)
                % Compute the interior angles by the cosine rule
                length4e = obj.length4e;
                weight = permute((ones(3) / 2 - eye(3)), [3 1 2]);
                obj.angle4e = squeeze(acos( ...
                    sum(bsxfun(@times, length4e.^2, weight), 2) ...
                    ./ prod(cat(3, length4e(:,[2 3]), length4e(:,[3 1]), ...
                                   length4e(:,[1 2])), 2)));
            end
            angle4e = obj.angle4e;
        end

        function angle = minAngle(obj)
            %%MINANGLE minimal interior angle of all triangles
            %(float: 1)
            angle = min(obj.angle4e(:));
        end

        function angle = maxAngle(obj)
            %%MAXANGLE minimal interior angle of all triangles
            %(float: 1)
            angle = max(obj.angle4e(:));
        end

        function length4ed = get.length4ed(obj)
            if isnan(obj.length4ed)
                % Compute the vector tangent to each edge
                U = obj.c4n(obj.n4ed(:,2),:) - obj.c4n(obj.n4ed(:,1),:);
                % Compute the length
                obj.length4ed = sqrt(sum(U.^2, 2));
            end
            length4ed = obj.length4ed;
        end

        function length4edDb = length4edDb(obj)
            %%LENGTH4EDDB length of the Dirichlet edges
            %(float: nDbEdges x 1)
            length4edDb = obj.length4ed(obj.DbEdges);
        end

        function length4edNb = length4edNb(obj)
            %%LENGTH4EDNB length of the Neumann edges
            %(float: nNbEdges x 1)
            length4edNb = obj.length4ed(obj.NbEdges);
        end

        function h4ed = get.h4ed(obj)
            if isnan(obj.h4ed)
                obj.h4ed = sqrt(obj.area4e(obj.e4ed(:,1)));
            end
            h4ed = obj.h4ed;
        end

        function length4e = length4e(obj, varargin)
            %%LENGTH4E length of the 3 edges of each triangle
            %(float: nElem x 3)
            %
            %   Note: This is necessary for Matlab's parallelisation.
            %         Rearrangement allows to address all edges for a given
            %         triangle.
            if nargin > 1
                length4e = obj.length4ed(obj.ed4e(varargin{:}));
            else
                length4e = obj.length4ed(obj.ed4e);
            end
        end

        function diam4e = diam4e(obj, varargin)
            %DIAM4E diameter of each triangle
            %(float: nElem x 1)
            diam4e = max(obj.length4e, [], 2);
        end

        function height4e = height4e(obj, varargin)
            %HEIGHT4E height of the triangle with respect to 
            %the three edges of each triangle
            %(float: nElem x 1)
            height4e = bsxfun(@rdivide, 2 * obj.area4e, obj.length4e);
        end

        %% GET FUNCTIONS FOR NODAL INFORMATION
        function c4e = c4e(obj)
            %%C4E coordinates of the three nodes on each element
            %(float: 2 x 3 x nElem)
            %
            %   Note: This is necessary for Matlab's parallelisation.
            %         Rearrangement allows to address all coordinates
            %         for a given triangle.
            c4e = permute(reshape(obj.c4n(obj.n4e,:), obj.nElem, 3, 2), ...
                          [3 2 1]);
        end

        %% GET FUNCTIONS FOR NEIGHBOURING INFORMATION
        function neighbours4e = get.neighbours4e(obj)
            % The order follows the local numbering of the edges
            % Contains 0 if it is a boundary edges
            if isnan(obj.neighbours4e)
                IndexSum4e = obj.e4ed(obj.ed4e,1) + obj.e4ed(obj.ed4e,2);
                obj.neighbours4e = reshape(IndexSum4e, [], 3)...
                                   - repmat((1:obj.nElem)', 1, 3);
            end
            neighbours4e = obj.neighbours4e;
        end

        %% GET FUNCTIONS FOR AFFINE TRANSFORMATIONS
        function Trafo4e = get.Trafo4e(obj)
            % This gives the Jacobian of the affine transformation from
            % the reference triangle to each triangle of the triangulation
            if isnan(obj.Trafo4e)
                obj.computeAffineTransformation();
            end
            Trafo4e = obj.Trafo4e;
        end

        function InvTrafo4e = get.InvTrafo4e(obj)
            % Inverse of Trafo4e
            if isnan(obj.InvTrafo4e)
                obj.computeAffineTransformation();
            end
            InvTrafo4e = obj.InvTrafo4e;
        end

        function Trafo4ed = get.Trafo4ed(obj)
            % This is Jacobian of the affine transformation from the
            % reference interval to each edge of the triangulation
            if isnan(obj.Trafo4ed)
                Trafo = ...
                    (obj.c4n(obj.n4ed(:,1),:) - obj.c4n(obj.n4ed(:,2),:))';
                obj.Trafo4ed = reshape(Trafo, [2 1 obj.nEdges]);
            end
            Trafo4ed = obj.Trafo4ed;
        end

        function Trafo4edDb = Trafo4edDb(obj)
            %%TRAFO4EDDB affine transformation for each Dirichlet edge
            %(int: 2 x 1 x nDbEdges)
            Trafo4edDb = obj.Trafo4ed(:,:,obj.DbEdges);
        end

        function Trafo4edNb = Trafo4edNb(obj)
            %%TRAFO4EDNB affine transformation for each Neumann edge
            %(int: 2 x 1 x nNbEdges)
            Trafo4edNb = obj.Trafo4ed(:,:,obj.NbEdges);
        end

        %% TRANSFORMATION FUNCTIONS
        function points4e = transform(obj, points, elements)
            %%TRANSFORM transforms the given points from the reference triangle
            %to each triangle given by the list elements
            %   MESHDATA.TRANSFORM(points, elements)
            %   MESHDATA.TRANSFORM(points) transforms points to each triangle of
            %       the triangulation
            if nargin < 3
                elements = 1:obj.nElem;
            end
            points = permute(points, [3 2 4 1]);
            points4e = sum(bsxfun(@times, obj.Trafo4e(:,:,elements), points), 2);
            position4e = permute(obj.c4n(obj.n4e(elements,3),:), [2 3 1]);
            points4e = bsxfun(@plus, points4e, position4e);
            points4e = permute(points4e, [1 3 4 2]);
        end

        function points = inverseTransform(obj, points4e, elements)
            %%INVERSETRANSFORM transforms the points from the triangles given
            %by the list elements to the reference triangle
            %   MESHDATA.TRANSFORM(points4e, elements)
            %   MESHDATA.TRANSFORM(points) transforms points from each triangle of
            %       the triangulation
            if nargin < 3
                elements = 1:obj.nElem;
            end
            position4e = permute(obj.c4n(obj.n4e(elements,3),:), [2 1]);
            points = permute(bsxfun(@minus, points4e, position4e), [4 1 2 3]);
            points = sum(bsxfun(@times, obj.InvTrafo4e(:,:,elements), points), 2);
            points = permute(points, [1 3 4 2]);
        end

        function points4ed = edgeTransform(obj, points, edges)
            %%EDGETRANSFORM transforms the given points from the reference edge
            %to each edge given by the list edges
            %   MESHDATA.EDGETRANSFORM(points, edges)
            %   MESHDATA.EDGETRANSFORM(points) transforms points to each
            %       edge of the triangulation
            if nargin < 3
                edges = 1:obj.nEdges;
            end
            points = permute(points, [3 2 4 1]);
            points4ed = sum(bsxfun(@times, obj.Trafo4ed(:,:,edges), points), 2);
            position4ed = permute(obj.c4n(obj.n4ed(edges,2),:), [2 3 1]);
            points4ed = bsxfun(@plus, points4ed, position4ed);
            points4ed = permute(points4ed, [1 3 4 2]);
        end

        %% GET FUNCTIONS FOR MID POINTS
        function mid4e = get.mid4e(obj)
            if isnan(obj.mid4e)
                obj.mid4e = (obj.c4n(obj.n4e(:,1),:)...
                             + obj.c4n(obj.n4e(:,2),:)...
                             + obj.c4n(obj.n4e(:,3),:)) / 3;
            end
            mid4e = obj.mid4e;
        end

        function mid4ed = get.mid4ed(obj)
            if isnan(obj.mid4ed)
                obj.mid4ed = (obj.c4n(obj.n4ed(:,1),:)...
                              + obj.c4n(obj.n4ed(:,2),:)) / 2;
            end
            mid4ed = obj.mid4ed;
        end

        function edgemid4e = edgemid4e(obj)
            %%EDGEMID4E midpoints of the three edges of each triangle
            %(float: 2 x 3 x nElem)
            edgemid4e = permute(reshape(obj.mid4ed(obj.ed4e,:), ...
                                        [obj.nElem 3 2]), [3 2 1]);
        end

        function mid4edDb = mid4edDb(obj)
            %%MID4EDDB midpoint of each Dirichlet edge
            %(float: nDbEdges x 2)
            mid4edDb = obj.mid4ed(obj.DbEdges,:);
        end

        function mid4edNb = mid4edNb(obj)
            %%MID4EDNB midpoint of each Neumann edge
            %(float: nNbEdges x 2)
            mid4edNb = obj.mid4ed(obj.NbEdges,:);
        end

        %% GET FUNCTIONS FOR NORMAL AND TANGENTIAL VECTOR FIELDS
        function normal4e = get.normal4e(obj)
            if isnan(obj.normal4e)
                obj.normal4e = ...
                    [obj.tangent4e(2,:,:); -obj.tangent4e(1,:,:)];
            end
            normal4e = obj.normal4e;
        end

        function normal4eDb = normal4eDb(obj)
            %%NORMAL4EDB outward unit normal vector of each Dirichlet edge
            %(float: nDbEdges x 2)
            normal4eDb = obj.normal4ed(obj.DbEdges,:) ...
                         .* repmat(obj.sign4edDb, 1, 2);
        end

        function normal4eNb = normal4eNb(obj)
            %%NORMAL4ENB outward unit normal vector of each Neumann edge
            %(float: nNbEdges x 2)
            normal4eNb = obj.normal4ed(obj.NbEdges,:) ...
                            .* repmat(obj.sign4edNb, 1, 2);
        end

        function tangent4e = get.tangent4e(obj)
            if isnan(obj.tangent4e)
                c4e = obj.c4e;
                obj.tangent4e = ...
                    bsxfun(@rdivide, c4e(:,[3 1 2],:) - c4e(:,[2 3 1],:),...
                           permute(obj.length4ed(obj.ed4e), [3 2 1]));
            end
            tangent4e = obj.tangent4e;
        end

        function normal4ed = get.normal4ed(obj)
            if isnan(obj.normal4ed)
                obj.normal4ed = [obj.tangent4ed(:,2), -obj.tangent4ed(:,1)];
            end
            normal4ed = obj.normal4ed;
        end

        function tangent4ed = get.tangent4ed(obj)
            if isnan(obj.tangent4ed)
                obj.tangent4ed = obj.c4n(obj.n4ed(:,2),:) ...
                                 - obj.c4n(obj.n4ed(:,1),:);
                obj.tangent4ed = bsxfun(@rdivide, obj.tangent4ed, obj.length4ed);
            end
            tangent4ed = obj.tangent4ed;
        end

        %% GET FUNCTIONS FOR SIGN OF NORMAL VECTORS
        function sign4e = get.sign4e(obj)
            % Contains the sign change of normal vectors of each
            % edge with the outward unit normal vector for every triangle
            if isnan(obj.sign4e)
                obj.computeEdges();
            end
            sign4e = obj.sign4e;
        end

        function sign4edDb = get.sign4edDb(obj)
            % Contains the sign change between edge normal and
            % outward normal of the domain
            if isnan(obj.sign4edDb)
                DbEdgeIndex = obj.ed4e(obj.e4edDb(:,1),:)...
                              == repmat(obj.DbEdges,1, 3);
                TMP = obj.sign4e(obj.e4edDb,:);
                obj.sign4edDb = TMP(DbEdgeIndex);
            end
            sign4edDb = obj.sign4edDb;
        end

        function sign4edNb = get.sign4edNb(obj)
            % Contains the sign change between edge normal and
            % outward normal of the domain
            if isnan(obj.sign4edNb)
                NbEdgeIndex = obj.ed4e(obj.e4edNb(:,1),:) ...
                              == repmat(obj.NbEdges,1, 3);
                TMP = obj.sign4e(obj.e4edNb,:);
                obj.sign4edNb = TMP(NbEdgeIndex);
            end
            sign4edNb = obj.sign4edNb;
        end

        %% GET FUNCTIONS FOR EDGE INFORMATION
        function n4ed = get.n4ed(obj)
            if isnan(obj.n4ed)
                obj.computeEdges();
            end
            n4ed = obj.n4ed;
        end

        function ed4e = get.ed4e(obj)
            if isnan(obj.ed4e)
                obj.computeEdges();
            end
            ed4e = obj.ed4e;
        end

        function e4ed = get.e4ed(obj)
            % Second entry is 0 for boundary edge
            if isnan(obj.e4ed)
                obj.computeEdges();
            end
            e4ed = obj.e4ed;
        end

        function edges = ed4n(obj, nodes)
            %%ED4N computes the edge numbers for every line in nodes
            %containing two nodes, returns 0 if the nodes form no edge
            %Due to memory reasons this is realised by a function
            %requiring a row-wise sorted version of n4ed
            %
            %   edges = ED4N(n4edSorted, nodes)
            if isnan(obj.n4edSorted)
                obj.n4edSorted = sort(obj.n4ed, 2);
            end
            [~, edges] = ismember(sort(nodes, 2), obj.n4edSorted, 'rows');
        end

        function isInteriorEdge = isInteriorEdge(obj)
            %%ISINTERIOREDGE provides boolean array determining whether a
            %given edge is an interior edge
            %(boolean: nEdges x 1)
            isInteriorEdge = (obj.e4ed(:,2) ~= 0);
        end

        function nInteriorEdges = nInteriorEdges(obj)
            %%NINTERIOREDGES returns number of interior edges
            nInteriorEdges = nnz(obj.isInteriorEdge);
        end

        function interiorEdges = interiorEdges(obj)
            %%INTERIOREDGES returns indices of interior edges
            Edges = (1:obj.nEdges)';
            interiorEdges = Edges(obj.isInteriorEdge);
        end

        function oppositeNode4ed = get.oppositeNode4ed(obj)
            if isnan(obj.oppositeNode4ed)
                obj.computeEdges();
            end
            oppositeNode4ed = obj.oppositeNode4ed;
        end

        function oppositeEdges4ed = get.oppositeEdges4ed(obj)
            if isnan(obj.oppositeEdges4ed)
                obj.computeEdges();
            end
            oppositeEdges4ed = obj.oppositeEdges4ed;
        end

        %% GET FUNCTIONS FOR INFORMATION ON DIRICHLET BOUNDARY
        function DbNodes = get.DbNodes(obj)
            if isnan(obj.DbNodes)
                obj.DbNodes = find(obj.onDirichlet(obj.c4n));
            end
            DbNodes = obj.DbNodes;
        end

        function n4edDb = n4edDb(obj)
            %%N4EDDB nodes of each Dirichlet boundary edge
            %(int: nDbEdges x 2)
            n4edDb = obj.n4ed(obj.DbEdges,:);
        end

        function c4nDb = c4nDb(obj)
            %%C4NDB coordinates of the Dirichlet nodes
            %(float: nDbNodes x 2)
            c4nDb = obj.c4n(obj.DbNodes,:);
        end

        function DbEdges = get.DbEdges(obj)
            if isnan(obj.DbEdges)
                obj.DbEdges = find(obj.onDirichlet(obj.mid4ed));
            end
            DbEdges = obj.DbEdges;
        end

        function e4edDb = e4edDb(obj)
            %%E4EDDB number of adjacent triangle for each Dirichlet edge
            %(int: nDbEdges x 1)
            e4edDb = obj.e4ed(obj.DbEdges, 1);
        end

        %% NEUMANN BOUNDARY
        function NbNodes = get.NbNodes(obj)
            if isnan(obj.NbNodes)
                obj.NbNodes = find(obj.onNeumann(obj.c4n));
            end
            NbNodes = obj.NbNodes;
        end

        function n4edNb = n4edNb(obj)
            %%N4EDNB nodes of each Neumann boundary edge
            %(int: nNbEdges x 2)
            n4edNb = obj.n4ed(obj.NbEdges,:);
        end

        function c4nNb = c4nNb(obj)
            %%C4NNB coordinates of the Neumann nodes
            %(float: nNbNodes x 2)
            c4nNb = obj.c4n(obj.NbNodes,:);
        end

        function NbEdges = get.NbEdges(obj)
            if isnan(obj.NbEdges)
                obj.NbEdges = find(obj.onNeumann(obj.mid4ed));
            end
            NbEdges = obj.NbEdges;
        end

        function e4edNb = e4edNb(obj)
            %%E4EDNB number of adjacent triangle for each Neumann edge
            %(int: nNbEdges x 1)
            e4edNb = obj.e4ed(obj.NbEdges, 1);
        end

        %% FUNCTIONS FOR MESH REFINEMENT
        function n4edMarked = markBulk(obj, eta4e, theta)
            %%MARKBULK marks a set of edges according to bulk criterion

            % Mark elements by bulk criterion
            [eta4e, Ind] = sort(eta4e, 'descend');
            cumsumEta4e = cumsum(eta4e);
            nMarked = ...
                find(cumsumEta4e >= theta * cumsumEta4e(end), 1, 'first');
            eMarked = Ind(1:nMarked);

            % Select edges of marked elements
            n4edMarked = obj.markEdges4Refinement(eMarked);
        end

        function n4edMarked = markMaximum(obj, eta4e, theta)
            %%MARKMAXIMUM marks a set of edges according to maximum criterion

            % Mark elements by bulk criterion
            eMarked = find(eta4e >= theta * max(eta4e));

            % Select edges of marked elements
            n4edMarked = obj.markEdges4Refinement(eMarked);
        end

        function n4edMarked = markNodes(obj, coordinates)
            %%MARKNODES marks a set of edges of the triangles that contain
            %at least one of nodes given by the coordinates

            % Mark elements that contain the given nodes
            eMarked = false(obj.nElem, 1);
            for j = 1:size(coordinates, 1)
                node = find(all(near(obj.c4n, coordinates(j,:)), 2));
                eMarked = (eMarked | any(near(obj.n4e, node), 2));
            end

            % Determine indices of marked elements
            eMarked = find(eMarked);

            % Select edges of marked elements
            n4edMarked = obj.markEdges4Refinement(eMarked);
        end

        function newMeshdata = refineUniformRed(obj)
            %%REFINEUNIFORMRED refine every triangle using the Red-strategy
            %   newMeshdata = MESHDATA.REFINEUNIFORMRED()

            %% REFINEMENT
            % Add new nodes to c4n
            c4n_new = [obj.c4n; obj.mid4ed];

            % Core routine: Create new elements
            %      n3       Element [n1 n2 n3] (as of a row in n4e)
            %     /  \      has the sides [s1 s2 s3] (as given by a
            %   s2    s1    row in s4e), in the order depicted here.
            %  /       \    Each side corresponds to a new node ---
            % n1 --s3-- n2  this node is the midpoint of the side.
            newNodes4e = obj.ed4e + obj.nNodes;
            n4e_new = [obj.n4e(:,1) newNodes4e(:,3) newNodes4e(:,2);  % [n1 s3 s2]
                       newNodes4e(:,3) obj.n4e(:,2) newNodes4e(:,1);  % [s3 n2 s1]
                       newNodes4e(:,2) newNodes4e(:,1) obj.n4e(:,3);  % [s2 s1 n3]
                       newNodes4e];                               % [s1 s2 s3]

            %% RECOVERING STRUCTURE FOR OLD DATA
            % Store relationship between old and new nodes
            parents4e_new = [1:obj.nElem 1:obj.nElem...
                             1:obj.nElem 1:obj.nElem];

            childNo4e_new = 4 * ones(4*obj.nElem, 1);

            % child numbers
            %       o
            %      / \
            %     / 3 \
            %    o-----o
            %   / \ 4 / \
            %  / 1 \ / 2 \
            % o-----o-----o
            obj.child4e = reshape(1:4*obj.nElem, [obj.nElem, 4]);

            obj.nChildren4e = 4 * ones(obj.nElem, 1);

            %% CREATE NEW MESHDATA OBJECT
            newMeshdata = MeshData(c4n_new, n4e_new, ...
                                   obj.onDirichlet, obj.onNeumann, ...
                                   parents4e_new, childNo4e_new);
        end

        function newMeshdata = refineBi3GB(obj, n4edMarked, allowHangingNodes)
            %%REFINEBI3GB refine using the Bisec3-Green-Blue-strategy
            %   newMeshdata = MESHDATA.REFINEBI3GB(n4edMarked) refines a given
            %   mesh using the Bisec3-Green-Blue refinement.

            %% PROCEED INPUT
            if nargin < 3
                allowHangingNodes = false;
            end

            %% CLOSURE
            n4edMarked = obj.closure(n4edMarked);

            %% UPDATE C4N
            mid4edMarked = ...
                (obj.c4n(n4edMarked(:,1),:) + obj.c4n(n4edMarked(:,2),:)) / 2;
            c4n_new = [obj.c4n; mid4edMarked];

            %% COMPUTE EDGES TO REFINE
            if allowHangingNodes
                % new nodes numbered according to n4sMarked
                [edge1Marked4e, edge1NewNodeNo4e] = ...
                    ismember(obj.n4e(:,[2 3]), n4edMarked, 'rows');
                [edge2Marked4e, edge2NewNodeNo4e] = ...
                    ismember(obj.n4e(:,[3 1]), n4edMarked, 'rows');
                [edge3Marked4e, edge3NewNodeNo4e] = ...
                    ismember(obj.n4e(:,[1 2]), n4edMarked, 'rows');
                isMarked4e = [edge1Marked4e, edge2Marked4e, edge3Marked4e];
                newNodes4e = zeros(nElem, 3);
                newNodes4e(isMarked4e) = ...
                    obj.nNodes + nonzeros([edge1NewNodeNo4e;
                                           edge2NewNodeNo4e;
                                           edge3NewNodeNo4e]);
            else
                MarkedEdges = obj.ed4n(n4edMarked);

                % newNodes4ed(k) = j > 0 if edge k is marked with new node j
                %                = 0     otherwise
                newNodes4ed = zeros(obj.nEdges, 1);
                newNodes4ed(MarkedEdges) = ...
                    obj.nNodes + (1:length(MarkedEdges));

                % newNodes4e(k,m) = j > 0 if edge m (1,2 or 3) of (old) element
                %                         k is marked with new node j
                %                 = 0     if edge m of (old) element k is not
                %                         marked
                newNodes4e = newNodes4ed(obj.ed4e);
            end

            %% REFINEMENT ALGORITHM
            %        n3        Element [n1 n2 n3] (as of a row in n4e)
            %       /  \       has the edges/new nodes [s1 s2 s3] (as
            %     s2    s1     given by a row in newNodes4e), in the order
            %    /        \    depicted here.  Each sj can be a new
            %   n1 --s3-- n2   node, depending of the refinement type.

            % Get lists of elements for corresponding refinement rules
            % Not to be refined
            Ind0 = ~any(newNodes4e,2);
            % To be bisec3 refined (all edges are marked)
            IndBi3 = all(newNodes4e,2);
            % To be blue_right refined (reference and 1st edge is marked)
            IndB1 = and(all(newNodes4e(:,[1 3]),2),~newNodes4e(:,2));
            % To be blue_left refined (reference and 2nd edge is marked)
            IndB2 = and(all(newNodes4e(:,[2 3]),2),~newNodes4e(:,1));
            % To be green refined (only reference edge is marked)
            IndG = and(newNodes4e(:,3),~any(newNodes4e(:,[1 2]),2));

            % Remaining cases for irregular refinement
            % Irregular refinement of edge 1 (green)
            IndIrr1 = and(newNodes4e(:,1),~any(newNodes4e(:,[2 3]),2));
            % Irregular refinement of edge 2 (green)
            IndIrr2 = and(newNodes4e(:,2),~any(newNodes4e(:,[1 3]),2));
            % Irregular refinement of edges 1 and 2 (blue right)
            IndIrr12 = and(all(newNodes4e(:,[1 2]),2),~newNodes4e(:,3));

            if any([IndIrr1; IndIrr2; IndIrr12])
                warning('Irregular refinement in refineBi3GB')
            end

            n4e_new = [% Untouched
                       obj.n4e(Ind0,:)
                       % Bisec3
                       %       o
                       %      /:\
                       %     / : \
                       %    o 2:3 o
                       %   / \ : / \
                       %  / 1 \:/ 4 \
                       % o-----o-----o
                       obj.n4e(IndBi3,1)     newNodes4e(IndBi3,3)  newNodes4e(IndBi3,2)
                       newNodes4e(IndBi3,3)  obj.n4e(IndBi3,3)     newNodes4e(IndBi3,2)
                       obj.n4e(IndBi3,3)     newNodes4e(IndBi3,3)  newNodes4e(IndBi3,1)
                       newNodes4e(IndBi3,3)  obj.n4e(IndBi3,2)     newNodes4e(IndBi3,1)
                       % Blue Right
                       %       o
                       %      /:\
                       %     / : \
                       %    /  :3 o
                       %   / 1 : / \
                       %  /    :/ 2 \
                       % o-----o-----o
                       obj.n4e(IndB1,3)    obj.n4e(IndB1,1)    newNodes4e(IndB1,3)
                       newNodes4e(IndB1,3) obj.n4e(IndB1,2)    newNodes4e(IndB1,1)
                       obj.n4e(IndB1,3)    newNodes4e(IndB1,3) newNodes4e(IndB1,1)
                       % Blue Left
                       %       o
                       %      /:\
                       %     / : \
                       %    o 2:  \
                       %   / \ : 3 \
                       %  / 1 \:    \
                       % o-----o-----o
                       obj.n4e(IndB2,1)    newNodes4e(IndB2,3) newNodes4e(IndB2,2)
                       newNodes4e(IndB2,3) obj.n4e(IndB2,3)    newNodes4e(IndB2,2)
                       obj.n4e(IndB2,2)    obj.n4e(IndB2,3)    newNodes4e(IndB2,3)
                       % Green
                       %       o
                       %      /:\
                       %     / : \
                       %    /  :  \
                       %   / 1 : 2 \
                       %  /    :    \
                       % o-----o-----o
                       obj.n4e(IndG,3) obj.n4e(IndG,1) newNodes4e(IndG,3)
                       obj.n4e(IndG,2) obj.n4e(IndG,3) newNodes4e(IndG,3)
                       % Irregular refinement of edge 1
                       %       o
                       %      / \
                       %     /   \
                       %    / 2  .o
                       %   /  .'   \
                       %  /.'   1   \
                       % o-----------o
                       obj.n4e(IndIrr1,1) obj.n4e(IndIrr1,2) newNodes4e(IndIrr1,1)
                       obj.n4e(IndIrr1,3) obj.n4e(IndIrr1,1) newNodes4e(IndIrr1,1)
                       % Irregular refinement of edge 2
                       %       o
                       %      / \
                       %     /   \
                       %    o.  2 \
                       %   /  ' .  \
                       %  /  1   ' .\
                       % o-----------o
                       obj.n4e(IndIrr2,1) obj.n4e(IndIrr2,2) newNodes4e(IndIrr2,2)
                       obj.n4e(IndIrr2,3) obj.n4e(IndIrr2,1) newNodes4e(IndIrr2,2)
                       % Irregular refinement of edge 1 and 2
                       %       o
                       %      / \
                       %     / 2 \
                       %    o----.o
                       %   / 3.'   \
                       %  /.'   1   \
                       % o-----------o
                       obj.n4e(IndIrr12,1)    obj.n4e(IndIrr12,2)    newNodes4e(IndIrr12,1)
                       newNodes4e(IndIrr12,1) obj.n4e(IndIrr12,3)    newNodes4e(IndIrr12,2)
                       obj.n4e(IndIrr12,1)    newNodes4e(IndIrr12,1) newNodes4e(IndIrr12,2)
                    ];

            %% RECOVERING STRUCTURE FOR OLD DATA
            Elements = permute(1:obj.nElem, [2 1]);
            nElem0 = nnz(Ind0);
            nElemBi3 = nnz(IndBi3);
            nElemB1 = nnz(IndB1);
            nElemB2 = nnz(IndB2);
            nElemG = nnz(IndG);
            nElemIrr1 = nnz(IndIrr1);
            nElemIrr2 = nnz(IndIrr2);
            nElemIrr12 = nnz(IndIrr12);

            parents4e_new = [repmat(Elements(Ind0),     1, 1)
                             repmat(Elements(IndBi3),   4, 1)
                             repmat(Elements(IndB1),    3, 1)
                             repmat(Elements(IndB2),    3, 1)
                             repmat(Elements(IndG),     2, 1)
                             repmat(Elements(IndIrr1),  2, 1)
                             repmat(Elements(IndIrr2),  2, 1)
                             repmat(Elements(IndIrr12), 3, 1)
                           ];

            childNo4e_new = [  ones(nElem0,     1)
                               ones(nElemBi3,   1)
                             2*ones(nElemBi3,   1)
                             3*ones(nElemBi3,   1)
                             4*ones(nElemBi3,   1)
                               ones(nElemB1,    1)
                             2*ones(nElemB1,    1)
                             3*ones(nElemB1,    1)
                               ones(nElemB2,    1)
                             2*ones(nElemB2,    1)
                             3*ones(nElemB2,    1)
                               ones(nElemG,     1)
                             2*ones(nElemG,     1)
                               ones(nElemIrr1,  1)
                             2*ones(nElemIrr1,  1)
                               ones(nElemIrr2,  1)
                             2*ones(nElemIrr2,  1)
                               ones(nElemIrr12, 1)
                             2*ones(nElemIrr12, 1)
                             3*ones(nElemIrr12, 1)
                            ];

            %% INFORMATION ON CHILDREN FOR CURRENT MESH
            obj.child4e = [ ...
                % unrefined
                repmat((1:nElem0)', 1, 4);...
                % bisec3
                (1:nElemBi3)'   + nElem0,...
                (1:nElemBi3)'   + nElem0 +   nElemBi3,...
                (1:nElemBi3)'   + nElem0 + 2*nElemBi3,...
                (1:nElemBi3)'   + nElem0 + 3*nElemBi3;...
                % blue_right
                (1:nElemB1)'    + nElem0 + 4*nElemBi3,...
                (1:nElemB1)'    + nElem0 + 4*nElemBi3,...
                (1:nElemB1)'    + nElem0 + 4*nElemBi3 +   nElemB1,...
                (1:nElemB1)'    + nElem0 + 4*nElemBi3 + 2*nElemB1;...
                % blue_left
                (1:nElemB2)'    + nElem0 + 4*nElemBi3 + 3*nElemB1,...
                (1:nElemB2)'    + nElem0 + 4*nElemBi3 + 3*nElemB1 +   nElemB2,...
                (1:nElemB2)'    + nElem0 + 4*nElemBi3 + 3*nElemB1 + 2*nElemB2,...
                (1:nElemB2)'    + nElem0 + 4*nElemBi3 + 3*nElemB1 + 2*nElemB2;...
                % green
                (1:nElemG)'     + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2,...
                (1:nElemG)'     + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2,...
                (1:nElemG)'     + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                                +   nElemG,...
                (1:nElemG)'     + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                                +   nElemG;...
                % irregular refinement of edge 1
                (1:nElemIrr1)'  + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                                + 2*nElemG,...
                (1:nElemIrr1)'  + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                                + 2*nElemG,...
                (1:nElemIrr1)'  + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                                + 2*nElemG +   nElemIrr1,...
                (1:nElemIrr1)'  + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                                + 2*nElemG +   nElemIrr1;...
                % irregular refinement of edge 2
                (1:nElemIrr2)'  + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                                + 2*nElemG + 2*nElemIrr1,...
                (1:nElemIrr2)'  + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                                + 2*nElemG + 2*nElemIrr1,...
                (1:nElemIrr2)'  + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                                + 2*nElemG + 2*nElemIrr1 +   nElemIrr2,...
                (1:nElemIrr2)'  + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                                + 2*nElemG + 2*nElemIrr1 +   nElemIrr2;...
                % irregular refinement of edge 1 and 2
                (1:nElemIrr12)' + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                                + 2*nElemG + 2*nElemIrr1 + 2*nElemIrr2,...
                (1:nElemIrr12)' + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                                + 2*nElemG + 2*nElemIrr1 + 2*nElemIrr2,...
                (1:nElemIrr12)' + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                                + 2*nElemG + 2*nElemIrr1 + 2*nElemIrr2 + nElemIrr12,...
                (1:nElemIrr12)' + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                                + 2*nElemG + 2*nElemIrr1 + 2*nElemIrr2 + nElemIrr12
                ];

            obj.nChildren4e = [  ones(nElem0,     1)
                               4*ones(nElemBi3,   1)
                               3*ones(nElemB1,    1)
                               3*ones(nElemB2,    1)
                               2*ones(nElemG,     1)
                               2*ones(nElemIrr1,  1)
                               2*ones(nElemIrr2,  1)
                               3*ones(nElemIrr12, 1)
                              ];

            %% CREATE NEW MESHDATA OBJECT
            newMeshdata = MeshData(c4n_new, n4e_new, ...
                                   obj.onDirichlet, obj.onNeumann, ...
                                   parents4e_new, childNo4e_new);
        end
    end

    %% PRIVATE AUXILIARY FUNCTIONS
    methods (Access = private)
        function obj = computeEdges(obj)
            %%COMPUTEEDGES computes all the edge information

            %% NODES FOR EDGES
            % Gather a list of all edges including duplicates
            % (occurring for inner edges), then sort each row and make
            % sure the rows are unique, thus eliminating duplicates.
            % The opposite nodes have the numbers 1, 2, 3. This
            % enables the common convention that node with number j
            % opposites the edge j.
            %
            %        2                      o
            %        | `                    | `
            %        |   `                  |   `
            %        |     `                |     `
            %        |       `              |       `
            %        |         `        -> [1]       [3]
            %        |           `          |           `
            %        |             `        |             `
            %        |               `      |               `
            %        3________________1     o______[2]_______o
            %

            % Create list of all edges
            [sortedEdges, localPosition] = sort([obj.n4e(:,[2 3]); ...
                                                 obj.n4e(:,[3 1]); ...
                                                 obj.n4e(:,[1 2])], 2);
            % Elininate duplicated edges
            %   NB: The following code ensures the convention that the
            %   nodes are ordered counter-clockwise with respect to T+
            [~, index, position] = unique(sortedEdges, 'rows', 'first');
            obj.n4ed = sortedEdges(sub2ind(size(sortedEdges), repmat(index, 1, 2), ...
                                      localPosition(index,:)));
            nEdges = size(index, 1);

            %% EDGES FOR ELEMENTS
            obj.ed4e = reshape(position, obj.nElem, 3);

            %% ELEMENTS FOR EDGES
            eNumbers = repmat((1:obj.nElem)', 3, 1);
            obj.e4ed = [eNumbers(index), ...
                        accumarray(position, eNumbers) - eNumbers(index)];

            %% INFORMATION ON INTERIOR EDGES
            isInteriorEdge = (obj.e4ed(:,2) ~= 0);
            Edges = (1:obj.nEdges)';
            interiorEdges = Edges(isInteriorEdge);
            nInteriorEdges = nnz(isInteriorEdge);

            %% NODES FOR EDGE PATCH
            obj.n4edgePatch = ...
                cat(3, obj.n4e(obj.e4ed(:,1),:), zeros(nEdges, 3, 1));
            obj.n4edgePatch(isInteriorEdge,:,2) = ...
                obj.n4e(obj.e4ed(isInteriorEdge,2),:);

            %% EDGE NUMBERS FOR EDGE PATCH
            obj.ed4edgePatch = ...
                cat(3, obj.ed4e(obj.e4ed(:,1),:), zeros(nEdges, 3, 1));
            obj.ed4edgePatch(isInteriorEdge,:,2) = ...
                obj.ed4e(obj.e4ed(isInteriorEdge,2),:);

            %% OPPOSITE NODES FOR EDGES
            % Local index of the opposite node in {1,2,3}
            LocalIndex = zeros(nEdges, 2);
            [LocalIndex(:,1),~] = ...
                find(all(bsxfun(@ne, obj.n4edgePatch(:,:,1), ...
                                permute(obj.n4ed, [1 3 2])), 3)');
            [LocalIndex(isInteriorEdge,2),~] = ...
                find(all(bsxfun(@ne, ...
                     obj.n4edgePatch(isInteriorEdge,:,2), ...
                     permute(obj.n4ed(isInteriorEdge,:), [1 3 2])), 3)');
            % Numbers of opposite nodes to each edge in T+ and T-
            obj.oppositeNode4ed = ...
                [obj.n4edgePatch(sub2ind(size(obj.n4edgePatch), ...
                                         Edges, LocalIndex(:,1))), ...
                 zeros(nEdges, 1)];
            obj.oppositeNode4ed(isInteriorEdge,2) = ...
                obj.n4edgePatch(sub2ind(size(obj.n4edgePatch), ...
                                        interiorEdges, ...
                                        LocalIndex(isInteriorEdge,2), ...
                                        2 * ones(nInteriorEdges, 1)));

            %% OPPOSITE EDGES FOR EDGES IN T+ AND T-
            SELECT = [2 3; 3 1; 1 2];
            obj.oppositeEdges4ed = ...
                cat(3, obj.ed4edgePatch(sub2ind(size(obj.ed4edgePatch), ...
                                          repmat(Edges, [1 2]), ...
                                          SELECT(LocalIndex(:,1),:))), ...
                    zeros(nEdges, 2));
            obj.oppositeEdges4ed(isInteriorEdge,:,2) = ...
                obj.ed4edgePatch(sub2ind(size(obj.ed4edgePatch), ...
                               repmat(interiorEdges, [1 2]), ...
                               SELECT(LocalIndex(isInteriorEdge,2),:), ...
                               2 * ones(nInteriorEdges, 2)));

            %% SIGN CHANGE OF NORMAL VECTORS
            obj.sign4e = ones(obj.nElem, 3);
            obj.sign4e(sub2ind(size(obj.sign4e), ...
                               obj.e4ed(isInteriorEdge,2), ...
                               LocalIndex(isInteriorEdge,2))) = -1;
        end

        function obj = computeAffineTransformation(obj)
            %%COMPUTEAFFINETRANSFORMATION computes Jacobian of inverse
            %affine transformation as well as its inverse

            % Initialise variables
            Trafo = zeros(2, 2, obj.nElem);
            InvTrafo = zeros(2, 2, obj.nElem);

            % Create local copies for use in parfor
            nElem = obj.nElem;
            c4e1 = permute(obj.c4n(obj.n4e(:,1),:), [2 1]);
            c4e2 = permute(obj.c4n(obj.n4e(:,2),:), [2 1]);
            c4e3 = permute(obj.c4n(obj.n4e(:,3),:), [2 1]);

            % Compute data
            parfor j = 1:nElem
                transformation = [c4e1(:,j) - c4e3(:,j),...
                                  c4e2(:,j) - c4e3(:,j)];
                Trafo(:,:,j) = transformation;
                % Inverting 2-by-2 matrix
                InvTrafo(:,:,j) = inv(transformation);
            end
            obj.Trafo4e = Trafo;
            obj.InvTrafo4e = InvTrafo;
        end

        function n4edMarked = closure(obj, n4edMarked)
            %%CLOSURE marks reference edges
            %   n4edMarked = MESHDATA.CLOSURE(n4edMarked) markes the reference
            %   edge of each element with at least one marked edge

            MarkedEdges = obj.ed4n(n4edMarked);
            isMarked = false(obj.nEdges, 1);
            isMarked(MarkedEdges) = true;
            isMarked4e = isMarked(obj.ed4e);

            %% CLOSURE
            while true
                % Check for consistency
                % isConsistent(k) == 0 if element k has a marked edge, but its
                %                      reference edge is not marked
                %                 == 1 otherwise
                isConsistent = ~and(any(isMarked4e(:,[1 2]),2), ~isMarked4e(:,3));
                if all(isConsistent)  % if all elements are consistent => done!
                    break;
                end
                % Mark edges that need to be refined
                n4edMarked_new = obj.n4e(~isConsistent,[1 2]);
                % Remove duplettes
                n4edMarked_new = unique(sort(n4edMarked_new, 2), 'rows', 'stable');
                n4edMarked = [n4edMarked; n4edMarked_new]; %#ok<AGROW>

                NewMarkedEdges = obj.ed4n(n4edMarked_new);
                isMarked(NewMarkedEdges) = 1;
                isMarked4e = isMarked(obj.ed4e);
            end
        end

        function n4edMarked = markEdges4Refinement(obj, eMarked)
            %%MARKEDGES4REFINEMENT marks all edges of the given elements for
            %refinement
            allEdgesMarked = ...
                [obj.n4e(eMarked,[1 2]);
                 obj.n4e(eMarked,[2 3]);
                 obj.n4e(eMarked,[3 1])];
            [~, IndMarked]   = unique(sort(allEdgesMarked, 2), 'rows');
            n4edMarked = allEdgesMarked(sort(IndMarked),:);
        end
    end

end

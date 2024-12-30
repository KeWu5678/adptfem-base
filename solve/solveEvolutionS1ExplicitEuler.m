function [u, ndof] = solveEvolutionS1ExplicitEuler(meshdata, A, b, r, f, I, u_0)
% A, b, r, f the data of the Parabolic problem 
% I the array of time step float x size(I)
% declare the global variables
nElem = meshdata.nElem;
nNodes = meshdata.nNodes;
area4e = meshdata.area4e;
c4e = meshdata.c4e;
n4e = meshdata.n4e; % indices of the 3 nodes of each triangle (int: nElem x 3)
c4n = meshdata.c4n;
mid4e = meshdata.mid4e;

% Given a coefficent vector x (float x nNodes) at time t, solve
% for the coeffient vector for dx\dt at time t.

% Abbreviate variables
n4e = meshdata.n4e;
area4e = meshdata.area4e;
c4e = meshdata.c4e;

% Evaluation of right-hand side f for midpoint quadrature rule,
% meshdata.mid4e(nElem x 2, barycenter)
f4e = f(meshdata.mid4e); % float x nElm
r4e = r(meshdata.mid4e); % float x nElm
b4e = b(meshdata.mid4e); % float x nElm x 2
A4e = A(meshdata.mid4e); % float x nElm x 2 x 2

% Pre-allocating variables
L4e = zeros(3, 3, nElem); % The mass matrix in front of the derivative
S4e = zeros(3, 3, nElem); % The local matrix that is multiplied by the vector
R4e = zeros(3, nElem); % RHS of the equation system for dx
cS4e = zeros(3, 3, nElem);

% Computation of local matrices
parfor j = 1:nElem
    grads = [ones(1,3); c4e(:,:,j)] \ [zeros(1,2); eye(2)];
    L4e(:,:,j) =  area4e(j) / 12 * (ones(3) + eye(3));
    S4e(:,:,j) = area4e(j) * (grads * ( A4e(:, :, j) * grads')) + ...
         grads * b4e(j, :)' * repmat(2/3 * area4e(j), 1, 3) + ...
        r4e(j) * area4e(j) / 12 * (ones(3) + eye(3));
    R4e(:,j) = f4e(j) * area4e(j) / 3 * ones(3, 1);
    cS4e(:,:,j) = area4e(j) * (grads * grads')
end

% Convert global into local degrees of freedom
ColumnIndices4e = permute(repmat(n4e, 1, 1, 3), [3 2 1]);
RowIndices4e = permute(ColumnIndices4e, [2 1 3]);
Indices4e = permute(meshdata.n4e, [2 1]);

% Assembling
LHS = sparse(RowIndices4e(:), ColumnIndices4e(:), L4e(:), ...
    nNodes, nNodes);
S = sparse(RowIndices4e(:), ColumnIndices4e(:), S4e(:), ...
    nNodes, nNodes);
RHS = accumarray(Indices4e(:), R4e(:), [nNodes 1]);

% Determining degrees of freedom
dof = setdiff(1:nNodes, meshdata.DbNodes);
ndof = length(dof);

    function dx = SolveEllip(x)
        % Solve linear system
        RHS_elli = RHS - S' * x;
        dx = zeros(nNodes, 1);
        dx(dof) = LHS(dof,dof) \ RHS_elli(dof);
    end


%% Naive explicit Euler method
u_ini = zeros(nNodes, 1);
for k = 1:nNodes
    u_ini(k) = u_0(c4n(k, 1), c4n(k, 2));
end

numSteps = size(I,2);
u = zeros(nNodes, numSteps);
du = zeros(nNodes, numSteps);
u(:, 1) = u_ini;
for k = 2: numSteps
    h = (I(k) - I(k - 1));
    dx = SolveEllip(u(:, k - 1));
    increment = h * dx;
    u_k = u(:, k -1) + increment;
    u(:, k) = u_k;
    du(:, k -1) = dx;
end
% debugging
% if isnan(err)
%     keyboard
% end
           
end


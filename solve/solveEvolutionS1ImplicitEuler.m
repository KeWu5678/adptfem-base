function [u, ndof] = solveEvolutionS1ImplicitEuler(meshdata, A, b, r, f, I, u_0)
% A, b, r, f the data of the Parabolic problem
% I the array of time step 1 x size(I)
% return u : nNodes x NumStepsize

%% Space discritisation
% declare the global variables
nElem = meshdata.nElem;
nNodes = meshdata.nNodes;
n4e = meshdata.n4e;
area4e = meshdata.area4e;
c4e = meshdata.c4e;
c4n = meshdata.c4n;
mid4e = meshdata.mid4e;

% Evaluation of right-hand side f for midpoint quadrature rule,
% meshdata.mid4e(nElem x 2, barycenter)
f4e = f(mid4e); % float x nElm
r4e = r(mid4e); % float x nElm
b4e = b(mid4e); % float x nElm x 2
A4e = A(mid4e); % float x nElm x 2 x 2

% Pre-allocating variables
M4e = zeros(3, 3, nElem); % The mass matrix in front of the derivative
S4e = zeros(3, 3, nElem); % The  matrix that is multiplied by the vector
R4e = zeros(3, nElem); % RHS of the equation system for dx

% Computation of local matrices
parfor j = 1:nElem
    grads = [ones(1,3); c4e(:,:,j)] \ [zeros(1,2); eye(2)];
    M4e(:,:,j) =  area4e(j) / 12 * (ones(3) + eye(3));
    S4e(:,:,j) = area4e(j) * (grads * ( A4e(:, :, j) * grads')) + ...
         grads * b4e(j, :)' * repmat(2/3 * area4e(j), 1, 3) + ...
        r4e(j) * area4e(j) / 12 * (ones(3) + eye(3));
    R4e(:,j) = f4e(j) * area4e(j) / 3 * ones(3, 1);
end

% Convert global into local degrees of freedom
ColumnIndices4e = permute(repmat(n4e, 1, 1, 3), [3 2 1]);
RowIndices4e = permute(ColumnIndices4e, [2 1 3]);
Indices4e = permute(meshdata.n4e, [2 1]);

% Assembling
M = sparse(RowIndices4e(:), ColumnIndices4e(:), M4e(:), ...
    nNodes, nNodes);
S = sparse(RowIndices4e(:), ColumnIndices4e(:), S4e(:), ...
    nNodes, nNodes);
R = accumarray(Indices4e(:), R4e(:), [nNodes 1]);

% Determining degrees of freedom
dof = setdiff(1:nNodes, meshdata.DbNodes);
ndof = length(dof);


%% Implement the Naive implicit Euler
% Basis transfomation of the initial condition
u_ini = zeros(nNodes, 1);
for k = 1:nNodes
    u_ini(k) = u_0(c4n(k, 1), c4n(k, 2));
end

numSteps = size(I,2);
u = zeros(nNodes, numSteps);
u(:, 1) = u_ini;
for k = 2: numSteps
    RHS = M * u(:, k - 1) + (I(k) - I(k - 1)) * R;
    MAT = M + (I(k) - I(k - 1)) * S';
    x = zeros(nNodes, 1);
    x(dof) = MAT(dof,dof) \ RHS(dof);
    u(:, k) = x; 
end
end


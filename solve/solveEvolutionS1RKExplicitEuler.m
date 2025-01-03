function [u, ndof] = solveEvolutionS1RKExplicitEuler(meshdata, A, b, r, f, I, u_0)
% A, b, r, f the data of the Parabolic problem 
% I the array of time step float x size(I)
% declare the global variables
nElem = meshdata.nElem;
nNodes = meshdata.nNodes;
c4n = meshdata.c4n;
mid4e = meshdata.mid4e;
n4e = meshdata.n4e;
area4e = meshdata.area4e;
c4e = meshdata.c4e;

%% Space discretisation

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

% Computation of local matrices
parfor j = 1:nElem
    grads = [ones(1,3); c4e(:,:,j)] \ [zeros(1,2); eye(2)];
    L4e(:,:,j) =  area4e(j) / 12 * (ones(3) + eye(3));
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
LHS = sparse(RowIndices4e(:), ColumnIndices4e(:), L4e(:), ...
    nNodes, nNodes);
S = sparse(RowIndices4e(:), ColumnIndices4e(:), S4e(:), ...
    nNodes, nNodes);
RHS = accumarray(Indices4e(:), R4e(:), [nNodes 1]);

% Determining degrees of freedom
dof = setdiff(1:nNodes, meshdata.DbNodes);
ndof = length(dof);

%% RK4 Implementation
    function [x_1] = RK4(x_0, h)
        k1 = -S * x_0;
        k2 = -S * (x_0 + 0.5 * h * k1);
        k3 = -S * (x_0 + 0.5 * h * k2);
        k4 = -S * (x_0 + h * k3);
        R_temp = LHS * x_0 + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*h;
        x_1(dof) = LHS(dof,dof) \ R_temp(dof);
    end

u_ini = zeros(nNodes, 1);
for k = 1:nNodes
    u_ini(k) = u_0(c4n(k, 1), c4n(k, 2));
end

numSteps = size(I,2);
u = zeros(nNodes, numSteps);
u(:, 1) = u_ini;

for k = 2: numSteps
    h = (I(k) - I(k - 1));
    u(:, k) = RK4(u(:, k -1), h);
end       
end


function [u, ndof, err_l2, err_h1] = solveEvolutionS1RKNImplicitEuler(meshdata, A, b, r, f, I, u_0, u_exact, RA, RB, RC)
% A, b, r, f the data of the Parabolic problem
% I the array of time step float x size(I)
% the stage of RK-method
% RA, RB, RC, coefficients of RK, 
% RA: s x s matrix, RB: s-dim vector RC: s-dim vector

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
    cS4e(:,:,j) = area4e(j) * (grads * grads')
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
cS = sparse(RowIndices4e(:), ColumnIndices4e(:), cS4e(:), ...
    nNodes, nNodes);
% Determining degrees of freedom
dof = setdiff(1:nNodes, meshdata.DbNodes);
ndof = length(dof);


%% The Implicit Runge-Kunta Method
% given u_0, h, calculate u_1
s = size(RA);
LHS_Dia = diag( 1 + diag(RA)); 
offDia_RA = RA - diag(diag(RA));
    function [u_1] = IRK(u_0, h)
        LHS = (M + h.* S') * (LHS_Dia)  + h.* S' * offDia_RA;
        RHS = - S' * u + R;
        k(dof) = LHS(dof) \ RHS(dof);
        u_1 = u_0 + RB .* k;
    end
 
% Basis transfomation of the initial condition
u_ini = zeros(nNodes, 1);
for l = 1:nNodes
    u_ini(l) = u_0(c4n(l, 1), c4n(l, 2));
end

numSteps = size(I,1);
u = zeros(nNodes, numSteps);
u(:, 1) = u_ini;

for l = 2: numSteps
    h = I(l) - I(l - 1);
    u(:, l) = RK2(u(:, l - 1), h, 0.001);
end

%% Error estimation 
% L2 error
err_l2 = 0;
for l = 2: numSteps
    for j = 1: nElem
        u_h = (u(n4e(j,1), l) + u(n4e(j,2), l) + u(n4e(j,3), l))/3;
        err_l2 = err_l2 + (u_exact(I(l),mid4e(j, 1), mid4e(j, 2)) - u_h)^2 * area4e(j) * (I(l) - I(l - 1));
    end
end
err_l2 = sqrt(err_l2);

%Calculate the error w.r.t. the H1 semi-norm
err_h1 = 0;
for l = 1: (numSteps - 1)
    % evaluation at each Nodes, reset every time step
    u_exact_h = zeros(nNodes, 1);
    for l = 1: nNodes
        u_exact_h(l) = u_exact(I(l), c4n(l, 1), c4n(l, 2)); 
    end
    err_h1 = err_h1 + (u_exact_h - u(:, l))' * S * (u_exact_h - u(:, l)) * (I(l + 1) - I(l));
end
err_h1 = sqrt(err_h1);


end


function [u, ndof] = solveEvolutionS1RKImplicitEuler(meshdata, A, b, r, f, I, u_0, Alpha, beta)
% A, b, r, f: the data of the Parabolic problem
% The coefficient of RK-method. 
% Alpha: s x s off-diagonal lower triangular matrix
% beta: s x 1 vector 
% I: the array of time step float 1 x numstep array

%% Space discritisation
% declare the global variables
nElem = meshdata.nElem;
nNodes = meshdata.nNodes; % nNodes x 1
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
L4e = zeros(3, 3, nElem); % The mass matrix in front of the derivative
S4e = zeros(3, 3, nElem); % The  matrix that is multiplied by the vector
R4e = zeros(3, nElem); % RHS of the equation system for dx

% Computation of local matrices
parfor j = 1:nElem
    grads = [ones(1,3); c4e(:,:,j)] \ [zeros(1,2); eye(2)];
    % In general the coefficient of dot{u} is the mass matrix
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


%% The Implicit Runge-Kunta Method
    function [x_1] = RK(x_0, h)
        % x_0 nNodes

        % Identify the stage
        s = size(beta, 1);

        % compute L_ij
        L = zeros(s * nNodes, s * nNodes);
        % the loop determine the submatrix L(j*s:j*s + nNodes, k*s: k*s +
        % nNodes)
        for j = 0:(s -1) 
            for k = 0: (s -1)
                if j +1 == k + 1
                    L((j*nNodes + 1):(j*nNodes + nNodes), (k*nNodes +1): (k*nNodes + nNodes)) = ...
                        (LHS + h .* S') .* (1 + Alpha(j + 1 , j + 1));
                else
                    L((j*nNodes + 1):(j*nNodes + nNodes), (k*nNodes +1): (k*nNodes + nNodes)) = ...
                         h * (sum(Alpha(j +1, :)) - Alpha(j + 1, j + 1)) * S';
                end
            end
        end

        % Compute k in the RK method(here denote by z)
        R = repmat(-S' * x_0 + RHS, s, 1);


        % Determining degrees of freedom
        % New set of all the nodes on the boundary
        dbnodes = meshdata.DbNodes;
        for j = 1:(s - 1)
            dbnodes_temp = meshdata.DbNodes + j * nNodes;
            dbnodes = cat(1, dbnodes, dbnodes_temp);
        end
        
        dof = setdiff(1:(nNodes * s), dbnodes); 
        ndof = length(dof);


        % Compute the k's of RK method. k1 ... ks stored as a vector
        z_temp = zeros(s * nNodes, 1);
        z_temp(dof) = L(dof, dof) \ R(dof);
        
        % Compute the x_1
        z = zeros(s, nNodes);
        for j = 1:s
            for l = 1:nNodes
                z(j, l) = z_temp((j - 1) * nNodes + l);
            end
        end
        x_1 = x_0 + h .* sum((beta .* z), 1)';
    end

%% Initialisation 

% Initialise the intial condition
u_ini = zeros(nNodes, 1);
for l = 1:nNodes
    u_ini(l) = u_0(c4n(l, 1), c4n(l, 2));
end

numSteps = size(I,2);
u = zeros(nNodes, numSteps);
u(:, 1) = u_ini;

for q = 2: numSteps
    h = (I(q) - I(q - 1));
    u(:, q) = RK(u(:, q -1), h);
end

ndof = ndof/ size(beta, 1);

end


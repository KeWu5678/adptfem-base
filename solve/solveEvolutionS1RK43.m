function [u, ndof, I] = solveEvolutionS1RK43(meshdata, A, b, r, f, T, u_0, tol, h_0)
% A, b, r, f the data of the Parabolic problem 
% T the terminal time

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

% Determining degrees of freedom
dof = setdiff(1:nNodes, meshdata.DbNodes); % set of all the nodes on the boundary 
ndof = length(dof);

%% The embedded Runge-Kutta Method(RK4(3)).
% Set up the Butcher tableau
Alpha = [0, 0, 0, 0, 0; ...
        1/2, 0, 0, 0, 0; ...
        0, 1/2, 0, 0, 0;...
        0, 0, 1, 0, 0;...
        1/6, 2/6, 2/6, 1/6, 0];
beta4 = [1/6; 2/6; 2/6; 1/6; 0];
beta3 = [1/6; 2/6; 2/6; 0; 1/6];

% Set up the parameter for the adaptive stepsize
alpha = 0.9;
h_min = 0.0001;
h_max = T;
a_max = 0.9;
a_min = 0.1;
s = 5;

% The embedded RK: Given a state x_0 and timestep h, calculate the next 
% state using 3rd order and 4th order RK. 
    function [x_1_3, x_1_4] = RK43(x_0, h)
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
        
        % Compute the RK4 and RK3
        z = zeros(s, nNodes);
        for j = 1:s
            for l = 1:nNodes
                z(j, l) = z_temp((j - 1) * nNodes + l);
            end
        end
        x_1_4 = x_0 + h .* sum((beta4 .* z), 1)';
        x_1_3 = x_0 + h .* sum((beta3 .* z), 1)';
    end

    % The error estimator
    function [x_1, t_1] = adpt_RK(x_0, t_0)
        [x_13, x_14] = RK43(x_0, h_0);
        eta = max(x_13 - x_14);
        %disp(eta);
        if eta <= tol
            t_1 = t_0 + h_0;
            x_1 = x_14;
        else
            q = min(a_max, max(a_min, alpha * (tol/eta)^(1/4)));
            h_1 = min([h_max, T - t_0, max([h_min, q * h_0])]);
            [~, x_1] = RK43(x_0, h_1);
            t_1 = t_0 + h_1;
        end
    end


%% The Iteration
% initialise the initial condition


u_ini = zeros(nNodes, 1);
for k = 1:nNodes
    u_ini(k) = u_0(c4n(k, 1), c4n(k, 2));
end

% final matrix to store the values
u = u_ini;
I = 0;

% initialise the iterator
t_0 = 0;
t_1 = 0;
u_0 = u_ini;

while t_0 < T
    [u_1, t_1] = adpt_RK(u_0, t_0);
    u = cat(2, u, u_1);
    I = cat(2, I, t_1);
    
    % iterate
    u_0 = u_1;
    t_0 = t_1;
end

end


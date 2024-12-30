function [err_l2, err_h1] = Error(meshdata, I, u, u_exact)
% The function calculate the L^infty(H) error of the Courant FEM method
% I: 1 x numSteps
% u: nNodes x numSteps
% u_exact: function handle

% declare the global variables
nElem = meshdata.nElem;
nNodes = meshdata.nNodes;
n4e = meshdata.n4e;
area4e = meshdata.area4e;
c4e = meshdata.c4e;
c4n = meshdata.c4n;
mid4e = meshdata.mid4e;

% Pre-allocating variables
S4e = zeros(3, 3, nElem);

% Computation of local matrices
parfor j = 1:nElem
    grads = [ones(1,3); c4e(:,:,j)] \ [zeros(1,2); eye(2)];
    S4e(:,:,j) = area4e(j) * (grads * grads')
end

% Convert global into local degrees of freedom
ColumnIndices4e = permute(repmat(n4e, 1, 1, 3), [3 2 1]);
RowIndices4e = permute(ColumnIndices4e, [2 1 3]);

% Assembling
S = sparse(RowIndices4e(:), ColumnIndices4e(:), S4e(:), ...
    nNodes, nNodes);

%% Calculate the L^infty(L^2) error
err_l2 = 0;
numSteps = size(I, 2);
for k = 2: numSteps
    % err_temp is the l2 error at time k using mitpoint evaluation
    for j = 1: nElem
        % 
        u_h = (u(n4e(j,1), k) + u(n4e(j,2), k) + u(n4e(j,3), k))/3;
        err_l2_temp = (u_exact(I(k), mid4e(j, 1), mid4e(j, 2)) - u_h)^2 * area4e(j);
    end
    % replace the err_l2 by err_l2_temp if applicable
    if err_l2 <= err_l2_temp
        err_l2 = err_l2_temp;
    end
end
err_l2 = sqrt( err_l2 * I(1, size(I, 2)));

%% Calculate the L^infty (H^1) error
err_h1 = 0;

for k = 2: numSteps
    % u_exact_temp is the evaluation at each Nodes, reset every time step
    u_exact_temp = zeros(nNodes, 1);
    for l = 1: nNodes
        u_exact_temp(l) = u_exact(I(k), c4n(l, 1), c4n(l, 2));
    end
    err_h1_temp = (u_exact_temp - u(:, k))' * S * (u_exact_temp - u(:, k));

    % set err_h1 = err_h1_temp if the latter is bigger
    if err_h1 <= err_h1_temp
        err_h1 = err_h1_temp;
    end
end

err_h1 = sqrt( err_h1* I(1, size(I, 2)));


    %% Initialisation
    normal4e = computeNormal4e(c4n,n4e);
    n4s = computeN4s(n4e);
    length4s = computeLength4s(c4n,n4s);
    nrElems = size(n4e,1);
    gradU = zeros(nrElems,2);
    jump4s = zeros(size(n4s,1),1);

    s4e = computeS4e(n4e);


    % %% Compute grad U
    % for elem = 1 : nrElems
    %     grads = [1,1,1;c4n(n4e(elem,:),:)']\[0,0;eye(2)];
    %     gradU(elem,:) = x(n4e(elem,:))' * grads;
    % end
    % 
    % %% inner jumps
    % for elem = 1 : size(n4e,1)
    %     jump4s(s4e(elem,:)) = jump4s(s4e(elem,:)) ...
    %                           + normal4e(:,:,elem)*gradU(elem,:)';
    % end
    % 
    % eta4s = (jump4s * length4s).^2;
end

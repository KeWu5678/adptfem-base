function main_adpt()
% This program solves the linear parabolic problems assuming the trivial
% initial contion and Dirichlet Boudary condition.

% Load mesh information
[c4n, n4e, onDirichlet, onNeumann] = Square();
meshdata = MeshData(c4n, n4e, onDirichlet, onNeumann);
 
% Defne right-hand side and righthand side
f = @(x) zeros(size(x, 1), 1);
r = @(x) zeros(size(x, 1), 1);
b = @(x) zeros(size(x, 1), 2);
A = @(x) repmat((1/pi) .* eye(2), 1, 1, size(x, 1));
u_0 = @(x_1, x_2) sin(pi * x_1) * sin(pi * x_2);
u_exact = @(t, x_1, x_2) exp(-2 * t) * sin(pi * x_1) * sin(pi * x_2);

k = 1;  % Iterator and bound of the ERR/NDOF set
N = 10; % num of iteration bound
ERR_L2 = zeros(N,1);
ERR_H1 = zeros(N,1);
num_steps = zeros(N, 1);
tol = zeros(N, 1);


%% Initialisation 
tic;
% Set up the meshdata
minNdof = 100;
while true
    nNodes = meshdata.nNodes;
    dof = setdiff(1:nNodes, meshdata.DbNodes);
    ndof = length(dof);
    % Break condition depending on number of degrees of freedom
    if ndof >= minNdof
        break
    end
    % Uniform refinement
    meshdata = meshdata.refineUniformRed();
end

% Iterate w.r.t. the tolerance
N = 15;
for k = 1:N
    % Solve Courant FEM
    tol_temp = 10^(-k * 5);
    [u, ndof, I] = solveEvolutionS1RK43(meshdata, A, b, r, f, 1, u_0, tol_temp, 0.1);   

    % Load the information
    [ERR_L2(k,:), ERR_H1(k, :)] = Error(meshdata, I, u, u_exact);
    num_steps(k) = size(I, 2);
    tol(k) = tol_temp;
    %display(I);
    printOutput(ndof, ERR_L2(k, :), ERR_H1(k, :), tol(k), size(I, 2));
end

% Plot
plotOutput(meshdata, u(:, end), tol, num_steps, 'num of steps');
toc;
end

%% Plot & Output functions
% The Plot function
function plotOutput(meshdata, u, tol, num_steps, name)
% Plot the function & error in one plot
% Choose figure
figure(1)

% Solution plot
subplot(1, 2, 1);
plotS1(meshdata, u);

subplot(1, 2, 2);
p = loglog(tol,num_steps);
p.LineStyle = "--";
p.Color = "red";
p.Marker = ".";
xlabel('tol')
ylabel('num of steps')
legend(name,'Location','northeast')
set(gca, 'xdir', 'reverse' )



% Finalise plots
drawnow;
end

% The print function
function printOutput(ndof, err1, err2, tol, I)
%%PRINTOUTPUT prints information to command line
fprintf('ndof= %d err_L2(L2)= %d  err_L2(H1)= %d tol = %d  size(I) = %d \n', ndof, err1,...
    err2, tol, I);
end

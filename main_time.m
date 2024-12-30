function main_time()
% This program solves the linear parabolic problems assuming the trivial
% initial contion and Dirichlet Boudary condition. 

% Load mesh information

[c4n, n4e, onDirichlet, onNeumann] = Square();
meshdata = MeshData(c4n, n4e, onDirichlet, onNeumann);

% Define right-hand side and righthand side
f = @(x) zeros(size(x, 1), 1);
r = @(x) zeros(size(x, 1), 1);
b = @(x) zeros(size(x, 1), 2);
A = @(x) repmat((1/(pi^2)) .* eye(2), 1, 1, size(x, 1));
u_0 = @(x_1, x_2) sin(pi * x_1) .* sin(pi * x_2);
u_exact = @(t, x_1, x_2) exp(-2 * t) .* sin(pi * x_1) .* sin(pi * x_2);

% The RK coefficients
Alpha4 = [0, 0, 0, 0;...
    0.5, 0, 0, 0;...
    0, 0.5, 0, 0;...
    0, 0, 1, 0];
beta4 = [1/6; 2/6; 2/6; 1/6];

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

% Start the Time Loop 
N = 10;
num_steps = zeros(N, 1);
ERR_L2 = zeros(N,1);
ERR_H1 = zeros(N,1);
for k = 1: N
    num_steps(k, :) = 2^k;
    I = linspace(0, 1, 2^k);

    % The Explicit Euler
    %         [u, ~] = solveEvolutionS1ExplicitEuler(meshdata,...
    %             A, b, r, f, I, u_0);

    % The Implicit Euler
    %         [u, ~] = solveEvolutionS1ImplicitEuler(meshdata,...
    %             A, b, r, f, I, u_0);

    % The Explicit RK4
    %         [u, ~] = solveEvolutionS1RKExplicitEuler(meshdata,...
    %             A, b, r, f, I, u_0);

    % The Implicit RK
    [u, ndof] = solveEvolutionS1RKImplicitEuler(meshdata, ...
        A, b, r, f, I, u_0, Alpha4, beta4);

    % Loard the information
    [ERR_L2(k, :), ERR_H1(k, :)] = Error(meshdata, I, u, u_exact);

    % Print the Information
    printOutput(num_steps(k, :), ERR_L2(k), ERR_H1(k));
end

plotOutput(meshdata, u(:, end), ERR_L2, ERR_H1, num_steps, 'C(L2)', 'C(H1)');
toc;
end

%% Plot & Output functions
function plotOutput(meshdata, u, err1, err2, numsteps, name_err1, name_err2)
    % Choose figure
    figure(1)

    % Solution plot
    subplot(1, 3, 1);
    plotS1(meshdata, u);

    subplot(1, 3, 2);
    loglog(numsteps, err1, "--")
    xlabel('num of steps')
    ylabel('error')
    legend(name_err1,'Location','northeast')

    subplot(1, 3, 3);
    loglog(numsteps, err2, "--")
    xlabel('num of steps')
    ylabel('error')
    legend(name_err2,'Location','northeast')

    % Finalise plots
    drawnow;   
end

function printOutput(h, errl2, errh1)
%%PRINTOUTPUT prints information to command line
    fprintf('numsteps = %d errl2 = %d errh1 = %d \n', h, errl2, errh1);
end




%% Spare codes
%     %Triangulation plot
%     subplot(1, 2, 1);
%     if meshdata.nElem < 1e4
%         plotTriangulation(meshdata);

%     loglog(num_steps, ERRH1(1,:,:),'--',num_steps, ERRH1(2,:,:),'-',num_steps, ERRH1(3,:,:), '--', num_steps, ERRH1(4,:,:), '-.')
%     xlabel('num of steps')
%     ylabel('error')
%     legend('ndof = 10', 'ndof = 100', 'ndof = 1000', 'ndof = 10000','Location','northeast')
%     drawnow; 
%     legend(name_err1,'Location','northeast')

    % Create a loop of time discretisation and observe the error/ndof relation
    
%     ERRL2 = zeros(5, N, 1);
%     ERRH1 = zeros(5, N, 1);
%     for i = 1:4
%         minNdof = 10^i;
%         while true
%             nNodes = meshdata.nNodes;
%             dof = setdiff(1:nNodes, meshdata.DbNodes);
%             ndof = length(dof);
%             % Break condition depending on number of degrees of freedom
%             if ndof >= minNdof
%                 break
%             end
%             % Uniform refinement
%             meshdata = meshdata.refineUniformRed();
%         end
%         err_l2 = zeros(N,1);
%         err_h1 = zeros(N,1);
%         num_steps = zeros(N, 1);
%         for k = 1: N
%             num_steps(k) = 2^k;
%             I = linspace(0, 1, 2^k);
% 
%             % Function & error value
%             [u, ~] = solveEvolutionS1RKExplicitEuler(meshdata,...
%                 A, b, r, f, I, u_0);
%             [err_l2(k), err_h1(k)] = Error(meshdata, I, u, u_exact);
% 
%             % Print the Information
%             printOutput(num_steps(k), err_l2(k), err_h1(k));
%         end
%         ERRL2(i, :, :) = err_l2;
%         ERRH1(i, :, :) = err_h1;
%     end

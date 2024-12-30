function main()
% Iteration of the meshsize 

    % Parameter for break condition(set by choice)
    minNdof = 1e3;

    % Load mesh information
    [c4n, n4e, onDirichlet, onNeumann] = Square();
    meshdata = MeshData(c4n, n4e, onDirichlet, onNeumann);

    % Define right-hand side and righthand side
    f = @(x) zeros(size(x, 1), 1);
    r = @(x) zeros(size(x, 1), 1);
    b = @(x) zeros(size(x, 1), 2);
    A = @(x) repmat((1/pi)^2 .* eye(2), 1, 1, size(x, 1));
    u_0 = @(x_1, x_2) sin(pi * x_1) * sin(pi * x_2);
    u_exact = @(t, x_1, x_2) exp(-2 * t) * sin(pi * x_1) * sin(pi * x_2);

    % The RK coefficients
    Alpha4 = [0, 0, 0, 0;...
            0.5, 0, 0, 0;...
            0, 0.5, 0, 0;...
            0, 0, 1, 0];
    beta4 = [1/6; 2/6; 2/6; 1/6];

    k = 1;  % Iterator and bound of the ERR/NDOF set
    N = 10; % num of iteration bound
    ERR_L2 = zeros(N,1);
    ERR_H1 = zeros(N,1);
    NDOF = zeros(N, 1);
    numsteps = 100;
    I = linspace(0, 1, numsteps);

   %% Start time measurement
    tic;
    while true
        % Explicit Euler
        % [u, ndof] =solveEvolutionS1ExplicitEuler(meshdata,...
        %   A, b, r, f, I, u_0);
        
        % the Implicit Euler
        [u, ndof] =solveEvolutionS1ImplicitEuler(meshdata,...
            A, b, r, f, I, u_0);

        % The Explicit RK4
%         [u, ndof] =solveEvolutionS1RKExplicitEuler(meshdata,...
%             A, b, r, f, I, u_0);

        % The Implicit RK
        % [u, ndof] = solveEvolutionS1RKImplicitEuler(meshdata, ...
        %     A, b, r, f, I, u_0, Alpha4, beta4);

        % Load the information
        [ERR_L2(k,:), ERR_H1(k, :)] = Error(meshdata, I, u, u_exact);
        NDOF(k) = ndof;
        printOutput(ndof, ERR_L2(k, :), ERR_H1(k, :));

        % Iterate
        k = k + 1;
        % Break condition depending on number of degrees of freedom
        if ndof >= minNdof || k >= N + 1
            break
        end

        % Uniform refinement
        meshdata = meshdata.refineUniformRed();
    end

    % End time measurement
    toc;
    
    % Plot
    plotOutput(meshdata, u(:, end), ERR_L2, ERR_H1, NDOF, 'C(L2)', 'C(H1)');
end


%% Plot & Output functions
% The Plot function
function plotOutput(meshdata, u, err1, err2, NDOF, name_err1, name_err2)
% Plot the function & error in one plot
    % Choose figure
    figure(1)

    % Solution plot
    subplot(1, 3, 1);
    plotS1(meshdata, u); 

    subplot(1, 3, 2);
    loglog(NDOF, err1, "--")
    xlabel('ndof')
    ylabel('error')
    legend(name_err1,'Location','northeast')

    subplot(1, 3, 3);
    loglog(NDOF, err2, "--")
    xlabel('ndof')
    ylabel('error')
    legend(name_err2,'Location','northeast')

    % Finalise plots
    drawnow; 
end

% The print function
function printOutput(ndof, err1, err2)
%%PRINTOUTPUT prints information to command line
    fprintf('ndof = %d ERRL2 = %d ERRH1 = %d \n', ndof, err1, err2);
end










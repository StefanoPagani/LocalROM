
% Main - local reduced-order model: convergence analysis 

%   Copyright (c) 2018, Politecnico di Milano 
%   localROM - Stefano Pagani <stefano.pagani at polimi.it>

clc
clear all
close all

% model parameters definition 
param(1) = 1;      % domain lenght
param(2) = 0.015;  % conducibility
param(3) = 0.5;    % recovery parameter
param(4) = 2;      % recovery parameter

% loop over the number of clusters N_c

for kclust = [ 4 6 8 10]
    
    %kclust

    % FN solver class constructor
    FNS = FNSolver(param, 1024, 0, 2, 400);

    % dimension of the testing set
    N_test = 25;

    % time step lenght
    dt = (FNS.tF-FNS.t0)/FNS.Nt;


    % loop over the POD tolerance

    tolvec = logspace(-1,-4,2);

    for itol = 1:2


        % LROM class constructor
        LROM = localizedReduction('PEBL',kclust);
        
        % offline procedure
        LROM = offline(LROM, 15, tolvec(itol));

        % for reproducibility
        rng('default')
        
        % loop over the test set
        for itest = 1:N_test

            % random realization
            ptest = 0.003 + rand(1,1)*(0.05-0.003);
            
            % uniform grid
            %ptest = 0.02 + (itest-1)/(N_test-1) *(0.05-0.02);

            % reduced solver
            [u,w] = FNS.solveROM(ptest, LROM);

            % full-order solver
            [uh,wh] = FNS.solveFOM(ptest);

            % \ell^2 error
            %err_u(itest) = dt * norm( u - uh, 2) ;
            
            % H1 error
            err_u(itest) = dt * sum( sqrt( diag( (u - uh)'*FNS.Xnorm*(u - uh) ) )./( 1 + sqrt( diag( (uh)'*FNS.Xnorm*(uh) ) ) ) ) ;


        end

        % mean error over the test set
        err_vec(itol) = mean(err_u);

        % computing the maximum number of basis functions among clusters
        n_basis(itol) = 0;
        for iC = 1:length(LROM.V)
            n_basis(itol) = max ( n_basis(itol) , size(LROM.V{iC},2) );
        end
    end

    % error plot
    semilogy( n_basis, err_vec, '-o', 'linewidth', 2 )
    hold all
    title('Reduction error local ROM')
    xlabel('$$max_{k=1,\ldots,N_c} n_k$$','Interpreter','Latex')
    ylabel('error')
    
    pause(1)
    
end

legend( '2', '4', '6', '8' , '10' )


clear all
close all
% construction of the reduced-order model 

for kclust = [2 4 6 8 10]
    
    kclust

   
%     LROM = offline(LROM, 50, 1e-8);
%     keyboard
    %LROM = localizedReduction('Time',kclust);

    % construction of the FOM model


    param(1) = 1; % domain lenght
    param(2) = 0.015;  % conducibility
    param(3) = 0.5;   % 
    param(4) = 2;

    FNS = FNSolver(param, 1024, 0, 2, 400);

    N_test = 50;

    dt = (FNS.tF-FNS.t0)/FNS.Nt;


    % loop over the POD tolerance
    tolvec = logspace(-1,-6,7);

    for itol = 1:7

        LROM = localizedReduction('kmeansParam',kclust);
        % offline procedure
        LROM = offline(LROM, 35, tolvec(itol));

        rng('default')
        % testing
        for itest = 1:N_test

            %ptest = 0.0035 + (itest-1)/(N_test-1) *(0.049-0.0035);
            ptest = 0.003 + rand(1,1)*(0.05-0.003);
            %ptest = 0.02 + (itest-1)/(N_test-1) *(0.05-0.02);

            
            [u,w] = FNS.solveROM(ptest, LROM);

            [uh,wh] = FNS.solveFOM(ptest);

            err_u(itest) = dt * norm( u - uh, 2) ;
            
%             if itest==28
%                 keyboard
%             end

        end

        err_vec(itol) = mean(err_u);

        n_basis(itol) = 0;
        for iC = 1:length(LROM.V)
            n_basis(itol) = max ( n_basis(itol) , size(LROM.V{iC},2) );
        end
    end


    semilogy( n_basis, err_vec, '-o', 'linewidth', 2 )
    hold all
    title('Reduction error')
    xlabel('n')
    ylabel('error')
    
end

legend( '2', '4', '6' , '8' , '10' )





% 
% % construction of the reduced-order model 
% 
% LROM = localizedReduction('kmeansState',4);
% 
% LROM = offline(LROM, 10, 1e-6);
% 
% 
% %%
% 
% param(1) = 1; % domain lenght
% param(2) = 0.015;  % conducibility
% param(3) = 0.5;   % 
% param(4) = 2;
% 
% FNS = FNSolver(param, 1024, 0, 2, 400)
% 
% [u,w] = FNS.solveROM(0.05, LROM);
% [uh,wh] = FNS.solveFOM(0.05);
% 
% norm(uh-u)
% 
% [X,Y] = meshgrid( linspace(0,FNS.L, FNS.Nh+1), linspace(FNS.t0,FNS.tF, FNS.Nt+1)  );
% 
% subplot(1,2,1)
% surf( X, Y , u')
% shading interp
% xlabel('x')
% ylabel('time')
% title('FOM solution: voltage')
% set(gca,'fontsize', 14)
% 
% subplot(1,2,2)
% surf( X, Y , w')
% shading interp
% xlabel('x')
% ylabel('time')
% title('FOM solution: recovery variable')
% set(gca,'fontsize', 14)
% 
% 

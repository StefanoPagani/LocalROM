clc
clear all


% construction of the reduced-order model 

LROM = localizedReduction('Global',4);

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

for jtol = 1:length(tolvec)

    for itol = 1:length(tolvec)

    
        % offline procedure
        LROM = offline(LROM, 35, tolvec(itol), tolvec(jtol));

        rng('default')
        % testing
        for itest = 1:N_test

            %ptest = 0.0035 + (itest-1)/(N_test-1) *(0.049-0.0035);
            ptest = 0.003 + rand(1,1)*(0.05-0.003);

            %ptest = 0.003 + (itest-1)/(N_test-1) *(0.05-0.003);

            [u,w] = FNS.solveROMHyperRed(ptest, LROM);

            [uh,wh] = FNS.solveFOM(ptest);

            err_u(itest) = dt * norm( u - uh, 2) ;

        end

        err_vec(itol) = mean(err_u);
        n_basis(itol) = size(LROM.V,2);

    end
    
    semilogy( n_basis, err_vec, '-o', 'linewidth', 2 )
    hold all
    title('Reduction error')
    xlabel('n')
    ylabel('error')
    
    pause(1)
    
end


legend('tol DEIM = 1e-1', 'tol DEIM = 1.47e-2', 'tol DEIM = 2.15e-3', 'tol DEIM = 3.16e-4' ...
    , 'tol DEIM = 4.64e-5' , 'tol DEIM = 6.81e-6' , 'tol DEIM = 1e-6')





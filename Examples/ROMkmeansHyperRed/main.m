clear all
close all
% construction of the reduced-order model 

for kclust = [5]
    
%     kclust

   
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
    
    tolvecDEIM = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6];

    for jtol = 1:length(tolvecDEIM)
        
        jtol

        for itol = 1:7

            LROM = localizedReduction('PEBL',kclust);
            % offline procedure
            LROM = offline(LROM, 35, tolvec(itol), tolvecDEIM(jtol));

%             keyboard
            
            rng('default')
            % testing
            for itest = 1:N_test

                %ptest = 0.0035 + (itest-1)/(N_test-1) *(0.049-0.0035);
                ptest = 0.003 + rand(1,1)*(0.05-0.003);
                %ptest = 0.02 + (itest-1)/(N_test-1) *(0.05-0.02);


                [u,w] = FNS.solveROMHyperRed(ptest, LROM);

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

        pause(1)
        
    end
end

legend('tol DEIM = 1e-2', 'tol DEIM = 1e-3', 'tol DEIM = 1e-4' , 'tol DEIM = 1e-5' , 'tol DEIM = 1e-6')

%legend( '2', '4', '6' , '8' , '10' )




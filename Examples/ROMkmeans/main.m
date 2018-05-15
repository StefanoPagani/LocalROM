% construction of the reduced-order model 

LROM = localizedReduction('Global',4);

LROM = offline(LROM, 10);


%%

param(1) = 1; % domain lenght
param(2) = 0.015;  % conducibility
param(3) = 0.5;   % 
param(4) = 2;

FNS = FNSolver(param, 1024, 0, 2, 400)

[u,w] = FNS.solveROM(0.05, LROM);

[X,Y] = meshgrid( linspace(0,FNS.L, FNS.Nh+1), linspace(FNS.t0,FNS.tF, FNS.Nt+1)  );

subplot(1,2,1)
surf( X, Y , u')
shading interp
xlabel('x')
ylabel('time')
title('FOM solution: voltage')
set(gca,'fontsize', 14)

subplot(1,2,2)
surf( X, Y , w')
shading interp
xlabel('x')
ylabel('time')
title('FOM solution: recovery variable')
set(gca,'fontsize', 14)



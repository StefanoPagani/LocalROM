% Main - Numerical approximation of the full-order model for an instance of the parameter epsilon

%   Copyright (c) 2018, Politecnico di Milano 
%   LocalROM - Stefano Pagani <stefano.pagani at polimi.it>

clc
clear all
close all

% model parameters definition
param(1) = 1;      % domain lenght
param(2) = 0.015;  % conducibility
param(3) = 0.5;    % recovery parameter
param(4) = 2;      % recovery parameter

% constructor
FNS = FNSolver(param, 1024, 0, 2, 400);

% solve the forward model
[u,w] = FNS.solveFOM(0.005);

% construct the grid for the figures
[X,Y] = meshgrid( linspace(0,FNS.L, FNS.Nh+1), linspace(FNS.t0,FNS.tF, FNS.Nt+1)  );

% figure 1: FOM voltage
%subplot(1,2,1)
figure()
surf( X, Y , u')
shading interp
xlabel('x')
ylabel('time')
title('FOM voltage ($$\epsilon$$=0.005)', 'Interpreter', 'LaTeX')
set(gca,'fontsize', 22)
axis([ 0 1 0 2 -0.5 1.5 ])

% figure 2: FOM recovery variable
%subplot(1,2,2)
figure()
surf( X, Y , w')
shading interp
xlabel('x')
ylabel('time')
%title('FOM solution: recovery variable (\epsilon=0.05)')
title('FOM recovery variable ($$\epsilon$$=0.005)', 'Interpreter', 'LaTeX')
set(gca,'fontsize', 22)
axis([ 0 1 0 2 0 0.2 ])


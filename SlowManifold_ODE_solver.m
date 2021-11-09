% Filename: SlowManifold_ODE_solver.m
% Author:  Maud El-Hachem
% Queensland University of Technology, Brisbane, Australia, November 2021
% Reference:  M. El-Hachem, S.W. McCue, M.J. Simpson (2021) 
% Travelling wave analysis of substrate-mediated cellular invasion
% This file script displays figure 8(b) that shows the vector field on the 
% slow manifold, superimosed with several solution trajectories.

clear;
% colours used
colors = [[0, 0.4470, 0.7410]; ...      
          	[0.8500, 0.3250, 0.0980;]];

% parameters r_1, r_2 and wave speed c
r1 = 1;
r2 = 1;
c = 10;

% displaying the uninvaded equilibrium point (0,0,0)
plot(0,0,'o','MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:),'LineWidth',4);
hold on

% domain for U and S comprised in the figure
Umin = 0;
Umax = 0.3;
Smin = 0;
Smax = 0.3;

% setting and computing the field vectors of the solution
U = linspace(Umin,Umax,10);
S = linspace(Smin,Smax,10);
[x,y] = meshgrid(U,S);
dU = zeros(size(x));
dS = zeros(size(x));
% ODE system for the vector field
f1 = @(X) [X(2).*(-X(1)/c + (c^2-r1)./c^3.*X(1).^2 -  (1-r2).*X(1).*X(2)./c^3);
    -r1/c*X(1).*X(2) + r2./c.*X(2).^2;];
for i = 1:numel(x)
    Xprime = f1([x(i); y(i)]);
    dU(i) = Xprime(1);
    dS(i) = Xprime(2);
end
% normalising the vector field
Un=dU./sqrt(dU.^2+dS.^2);
Sn=dS./sqrt(dU.^2+dS.^2);
% showing the vector field
q = quiver(x,y,Un,Sn,'k','LineWidth',.75);
set(q,'Autoscale','on', 'AutoscaleFactor',.5)

% computing and showing several trajectories solutions of the ODE system
Us=Umax;
for Ss=0.01:0.01:Smax
    % calling function heunsolver with wave speed c, parameters r1 and r2,
    % dzeta = 0.01, zeta going from 0 to 1000, and the initial conditions
    % U(0) = Us and S(0) = Ss
    [U, S] = heunSolver(c,r1,r2,0.01,0,1000,Us,Ss);
    % displaying the solution from the ode
    plot(U,S,'Color',colors(2,:),'LineWidth',1);
end

% displaying the labels of the corresponding axis
xlabel('$U(\zeta)$','interpreter','latex','fontsize',18);
ylabel('$S(\zeta)$','interpreter','latex','fontsize',18);

% formatting the figure
xlim([Umin,Umax]);
ylim([Smin,Smax]);
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(gca,'fontsize',18);
box on


% Function heunSolver
% This function solves Equations (26) and (28) by Heun's method 
% as described in the supplementary material.
% INPUT ARGUMENTS:
% ** c, the wave speed
% ** dzeta, the step size used to discretise the domain of zeta
% ** z_begin and z_end, the lower and upper limit of the numerical domain 
% of zeta such as z_begin <= z <= z_end. The initial conditions are applied at
% zeta = z_begin.
% ** Uinit, Sinit, the values of the initial conditions
% OUTPUT ARGUMENTS:
% ** U : The solution U(zeta)
% ** S : The solution S(zeta)
% The output size for arrays U and S may not be equal to the original number
% of nodes correponding to (z_end-z_begin)/dz+1 in the discretised domain.
function [U, S] = heunSolver(c,r1,r2,dzeta,z_begin,z_end,Uinit,Sinit)
    zeta = z_begin:dzeta:z_end;
    % number of nodes in the domain
    sz = length(zeta);

    % initialisation 
    U = zeros(sz,1);
    S = zeros(sz,1);

    U(1) = Uinit;
    S(1) = Sinit;

    szout = sz;
    
	% main loop of Heun's method
    for i=1:sz-1
        Ubar = U(i) + dzeta*S(i)*(-U(i)/c + (c^2-r1)/c^3*U(i)^2 -(1-r2)* U(i)*S(i)/c^3);
        Sbar = S(i) + dzeta*(-r1/c*U(i)*S(i) + r2/c*S(i)^2);
        U(i+1) = U(i) + dzeta/2*(S(i)*(-U(i)/c + (c^2-r1)/c^3*U(i)^2 -(1-r2)*  U(i)*S(i)/c^3) ...
        + Sbar*(-Ubar/c + (c^2-r1)/c^3*Ubar^2  -(1-r2)*Ubar*Sbar/c^3));
        S(i+1) = S(i) + dzeta/2*(-r1/c*U(i)*S(i) + r2/c*S(i)^2 + -r1/c*Ubar*Sbar + r2/c*Sbar^2);

		% if the solution is outside the domain 0<U<1 and 0<S<r1/r2
        if (U(i+1) < 0 || S(i+1) > 1 || U(i+1) < 0 || S(i+1) > r1/r2)
            szout = i;
            break;
        end
    end
  
    % output variable
    U = U(1:szout,1);
    S = S(1:szout,1);

end
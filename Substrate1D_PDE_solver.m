% Filename: Substrate1D_PDE_solver.m
% Authors: Matthew Simpson and Maud El-Hachem
% Queensland University of Technology, Brisbane, Australia, November 2021
% Reference:  M. El-Hachem, S.W. McCue, M.J. Simpson (2021) 
% Travelling wave analysis of substrate-mediated cellular invasion
% This file script contains:
%   - two calls to function PDE_solver that solves Eqs (6)-(8) using 
%   Newton-Raphson method with one of the two possible initial conditions:
%   compact support as in Eqs (9)-(10) or exponential decay front
%   as in Eqs (11)-(12). The function displays time-dependent
%   solutions u(x,t) and s(x,t) as in Figures 3(b) and (f).
%	- function PDE_solver
%   - function tridia

% Generating Figure 3(b)
% Calling function PDE_solver(r1,r2,decayrate,T,L,dx,dt,tprint)
% with the following parameters: r1 = 10, r2 = 20, decayrate set to Inf, 
% meaning that we are not using the IC condition with exponential decay, 
% T = 60, L = 80, dx = 0.01, dt = 0.001 and array of times at which 
% the solutions u(x,t) and s(x,t) are printed, [0 20 40 60]
PDE_solver(10, 20, Inf, 60, 80, 0.01, 0.001, [0 20 40 60]);

% Generating Figure 3(f)
% Calling function PDE_solver(r1,r2,decayrate,T,L,dx,dt,tprint)
% with the following parameters: r1 = 10, r2 = 20, decayrate = 0.5, such 
% as the expected final speed is 2, T = 60, L = 200, dx = 0.01, 
% dt = 0.001 and array of times at which the solutions u(x,t) and s(x,t)
% are printed, [0 20 40 60]
PDE_solver(10, 20, 0.5, 60, 200, 0.01, 0.001, [0 20 40 60]);

% Function PDE_solver
% This function solves Eqs (6)-(8) using Newton-Raphson method 
% with one of the two possible initial conditions:
% compact support as in Eqs (9)-(10) or exponential decay front
% as in Eqs (11)-(12). The function displays time-dependent
% solutions u(x,t) and s(x,t) at requested times.
% INPUT ARGUMENTS:
% -- r1 and r2: parameters r1 and r2 from Eq (7), such as r1 > 0 and r2 > 0
% -- decayrate: if the initial condition used include a front with an
% exponential decay, decayrate represent the exponential decay rate a as in
% Eq 11. If the initial condition used is on compact support, decayrate
% must be set to Inf
% -- L: total length of domain
% -- T: final time of solution, starting at t=0
% -- dx: spacing between two nodes in spatial grid
% -- dt: time step
% -- tprint: a vector of time steps at which the solutions are printed
% OUPUT ARGUMENT: ws: the wave speed of the moving front measured using a 
% contour value u* = 0.5. 
function ws = PDE_solver(r1,r2,decayrate,T,L,dx,dt,tprint)

    % colors used for printing the density profiles
    colors = [0 0.4470 0.7410; 0.4660 0.6740 0.1880;0.9290 0.6940 0.1250;];
    % minimum tolerance used in Newton-Raphson algorithm
    tol=1e-6;
    % total number of time steps
    maxsteps=T/dt;
    % total number of nodes
    N=L/dx+1;
    % array for density profile u(x,t) at time step i+1
    u=zeros(1,N);
    % parameter beta in initial conditions 
    beta = 10;
    for i=1:N
        xx(1,i)=0.0+dx*(i-1);
    end
    
    % intial conditions for u and s
    % compact support
    if isinf(decayrate)
        for i=1:N
            if (xx(1,i) < 10)
                u(1,i) = 1;
            end
        end
    else
        % exponential decay
        for i=1:N
            if (xx(1,i) < beta)
                u(1,i) = 1;
            else
                u(1,i) = 1*exp(-decayrate*(xx(1,i)-beta));
            end
        end
    end
    
    % array for density profile u(x,t) at time step i
    pu=u;
    % array for density profile s(x,t) at time step i+1
    s=zeros(1,N);
    % array for density profile s(x,t) at time step i
    ps=s;

    % array for recording the position of the moving front at u=0.5
    Lrecord=zeros(maxsteps+1,1);
    Lrecord(1,1)=0;

    % coefficients used in the tridiagonal matrix
    a=zeros(1,N);
    b=zeros(1,N);
    c=zeros(1,N);
    d=zeros(1,N);

    % initial time
    t=0.0;
    % showing the initial conditions in a new figure
    figure;
    hold on
    plot(xx,u,'LineWidth',2,'Color',colors(1,:));
    plot(xx,s,'LineWidth',2,'Color',colors(2,:));
    
    % Newton-Raphson solver
    for i=1:maxsteps
        t=t+dt;
        kk=0;
        delu=ones(1,N);
        while norm(delu,Inf) > tol
            kk=kk+1;
    
        % coefficients from bondaries conditions for u
        a(1,1)=0;
        b(1,1)=-1.0;
        c(1,1)=1.0;
        d(1,1)=-1*(u(1,2)-u(1,1));

        a(1,N)=0.0;
        b(1,N)=1.0;
        c(1,N)=0.0;
        d(1,N)=-1*(u(1,N));

        % coefficients from PDE du/dt
        for j=2:N-1
            a(1,j)=(s(1,j-1)+s(1,j))/(2*dx^2);
            b(1,j)=-1.0/dt-(s(1,j-1)+2*s(1,j)+s(1,j+1))/(2*dx^2)-(1-2*u(1,j));
            c(1,j)=(s(1,j+1)+s(1,j))/(2*dx^2);
            d(1,j)=(u(1,j)-pu(1,j))/dt-((s(1,j+1)+s(1,j))*(u(1,j+1)-u(1,j))-(s(1,j-1)+s(1,j))*(u(1,j)-u(1,j-1)))/(2*dx^2)-u(1,j)*(1-u(1,j));
        end

        delu = thomas(N,a,b,c,d);

        u(1,:)=u(1,:)+delu(1,:); 

        % coefficients from PDE ds/dt
        for j=1:N
            a(1,j)=0;
            b(1,j)=-1.0/dt;
            c(1,j)=0;
            d(1,j)=(s(1,j)-ps(1,j))/dt-r1*u(1,j)+r2*s(1,j);
        end
        deln = thomas(N,a,b,c,d);
        s(1,:)=s(1,:)+deln(1,:); 

        end

        % showing current solutions u and s at required time steps
        if (isempty(find(tprint==i*dt,1))==false)
            plot((0:N-1)*dx,u,'LineWidth',2,'Color',colors(1,:));
            plot((0:N-1)*dx,s,'LineWidth',2,'Color',colors(2,:));
            fprintf('Time %d\n',t); 
            fprintf('Iteration %d\n',kk);
        end

        % contour at u=0.5
        ii=find(u<=0.5,1);
        xii=(ii-1)*dx;
        xim=xii-1;

        % estimate of the current position of the front using a contour
        % at u=0.5
        xf=(0.5*dx+u(1,ii)*xim-u(1,ii-1)*xii)/(u(1,ii)-u(1,ii-1));

        % recording the current position of the front at u=0.5
        Lrecord(i+1,1)=xf;
        
        % updating the current arrays u and s
        pu(1,:)=u(1,:);
        ps(1,:)=s(1,:);
    end
    
    % formattig the figure
    ylabel('$u(x,t) \quad s(x,t)$','interpreter','latex');
    xlabel('$x$','interpreter','latex'); 
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(gca,'fontsize', 24);
    box on

    text(50,0.75,strcat(strcat('$r1=',num2str(r1,4)),'$'),'interpreter','latex','fontsize',24)
    text(50,.25,strcat(strcat('$r2=',num2str(r2,4)),'$'),'interpreter','latex','fontsize',24)
  
    % calculating the wave speed by interpolating the variabtion of
    % the recorded positions at u=0.5 for the time steps since t=T/2
    % we suppose that the speed has reached a constant level since t=T/2
    tt=0:dt:T;
    tt_end = tt(end-round(0.5*length(tt))+1:end);
    LL_end = Lrecord(end-round(0.5*length(tt))+1:end);
    hold on
    plot(tt_end,LL_end,'ro--')
    c=polyfit(tt_end,LL_end,1);
    format long
    % the wave speed
    ws = c(1);
end

%% Function tridia
% This function implements Thomas algorithm that solves a system
% of equations Ax = d, where A is a tridiagonal matrix. The parameters 
% a,b and c are the three diagonals of the matrix A. N is the size of 
% the vector solution x.
function x = thomas(N,a,b,c,d)
    x=zeros(1,N);
    bb=b;
    dd=d;
    for i=2:N
        ff=a(i)/bb(i-1);
        bb(i)=bb(i)-c(i-1)*ff;
        dd(i)=dd(i)-dd(i-1)*ff;
    end
    
    for i=1:N-1
    x(N)=dd(N)/bb(N);    
    j=N-i;
    x(j)=(dd(j)-c(j)*x(j+1))/bb(j);
    end
end




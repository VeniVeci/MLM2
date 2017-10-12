% Solver_Demo
% Burgers equation 1D case, finite element method, Dirichlet boundary
% conditions & homeogeneous B.C, solver.Function.

% Problem model:
% [u(t,x)]_t+1/2[u^2(t,x)]_x-q [u(t,x)]_xx=f(t,x), x \in (0,1),t>0;

% Boundary Condition (Dirichlet):
% u(t,0)=0; u(t,1)=0;

% Initial Conditions:
% u(0,x)=u_0(x);

% close all
clear;

%% Setup
% ------------------Problem Parameters------------------------------------- 
Paras.Re=1000;    % Reynolds Number
% v=1/Re;         % viscosity

% ------------------Solver Parameters--------------------------------------
Paras.n=128;           % Total Spatial elements
Paras.t_end=10;        % End time
Paras.t_n=500;        % Number of time step %200
% Paras.t=0:(t_end/t_n):t_end; % time sequence (t=0 does not account 1 step)

% solver = 'ode45';
% options = odeset('RelTol',1e-6,'AbsTol',1e-10);

% ------------------Calculating temp variable----------------------------- 
% h=1/n;      % space step size
% x = 0:h:1;  % coordinate sequence

%Initial condition
u0a=5;
u0b=1; 
u0=@(x) 1*exp(-(u0a*x+u0b)).*sin(pi*x);

%Sources term
g=@(x) 0.001*exp(x);


%% Main 
[Y,T,Time_Ode_solver]=burgerSolver(Paras,g,u0);


% Ploting
figure(1)
meshc(Y)
title(sprintf('Full FEM Solution'));


h=1/Paras.n;      % space step size
x = 0:h:1;  % coordinate sequence
x = 0:1/Paras.n:1;  % coordinate sequence
Y_Maxi=max(Y(:));
Y_Mini=min(Y(:));

for i = 1:Paras.t_n/100:Paras.t_n+1
    figure(3)
    plot(x,Y(:,i))
    axis([0,1,Y_Mini,Y_Maxi]);
    title(sprintf('Animation t= %0.3f',((i-1)*(Paras.t_end/Paras.t_n))))
    pause(0.1);  
end



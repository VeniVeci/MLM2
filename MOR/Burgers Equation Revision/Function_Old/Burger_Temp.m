%% Temperatury filed for burger solver 


%% Setup
% ------------------Problem Parameters------------------------------------- 
Paras.u0a=5; 
Paras.u0b=0.2; 

problemPara.re=1000; % Reynolds Number
% v=1/Re;         % viscosity

% ------------------Solver Parameters--------------------------------------
solverPara.n=128;           % Total Spatial elements
solverPara.t_end=10;        % End time
solverPara.t_n=500;        % Number of time step %200
% Paras.t=0:(t_end/t_n):t_end; % time sequence (t=0 does not account 1 step)

solverPara.solver = 'ode45';
% solverPara.options = odeset('RelTol',1e-6,'AbsTol',1e-10);

% ------------------Calculating temp variable----------------------------- 
% h=1/n;      % space step size
% x = 0:h:1;  % coordinate sequence

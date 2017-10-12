% burgerPod_deimSolver_Demo
% Burgers equation 1D case, finite element method, Dirichlet boundary
% conditions & homeogeneous B.C, solver.

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
% v=1/Re;     % viscosity

% ------------------Solver Parameters--------------------------------------
Paras.n=128;           % Total Spatial elements
Paras.t_end=4;        % End time
Paras.t_n=100;        % Number of time step 
% Paras.t=0:(t_end/t_n):t_end; % time sequence (t=0 does not account 1 step)

% solver = 'ode45';
% options = odeset('RelTol',1e-6,'AbsTol',1e-10);

approximate_degree=10;
approximate_degree_DEIM=5;


%Initial condition
u0a=5; 
u0b=1; 
u0=@(x) 1*exp(-(u0a*x+u0b)).*sin(3*pi*x);
 
%Sources term
gx=@(x) 0.02*exp(x);


%% Main 

% Normal solver
[Y1,T1,TimeCost_FOM]=burgerSolver(Paras,gx,u0);

% POD via SVD
Y1_inter=Y1(2:end-1,:);
[U,S,~]=svd(Y1_inter,'econ');  % U*S*V'=Rec_X
U=U(:,1:approximate_degree);
eigenvalues=diag(S);

% POD MOR
% [Y2,T2,Time_Ode_MORsolver]=Burger1D_FEM_DBC_MOR_SolverF(Paras,U);
[Y2,T2,TimeCost_POD]=burgerPodSolver(Paras,gx,u0,U);

% POD on non-linear terms. This is only valid for the Burgers equation
% where the nonlinear term is element-wise square of u(x)
[U_DEIM,S_DEIM,~]=svd(Y1_inter.^2); 
[~,U_DEIM,P] = DEIM(U_DEIM);

[Y3,T3,TimeCost_DEIM_POD]=burgerPod_deimSolver(Paras,gx,u0,U,U_DEIM(:,1:approximate_degree_DEIM),P(:,1:approximate_degree_DEIM));




% Ploting
figure(1)
meshc(Y1)
title(sprintf('FEM Solution. Time cost= %0.3f',TimeCost_FOM))

figure(2)
meshc(Y2)
title(sprintf('MOR FEM Solution. Time cost= %0.3f',TimeCost_POD))

figure(3)
meshc(Y3)
title(sprintf('MOR FEM Solution. Time cost= %0.3f',TimeCost_DEIM_POD))


h=1/Paras.n;      % space step size
x = 0:h:1;  % coordinate sequence
for i = 1:Paras.t_n+1
    figure(4)
    
    subplot(3,1,1)
    plot(x,Y1(:,i))
%     axis([0,1,0,1]);
    title(sprintf('FEM Animation t= %0.3f',((i-1)*(Paras.t_end/Paras.t_n))))
    
    subplot(3,1,2)
    plot(x,Y2(:,i))
%     axis([0,1,0,1]);
    title(sprintf('MOR FEM Animation t= %0.3f',((i-1)*(Paras.t_end/Paras.t_n))))   
    
    subplot(3,1,3)
    plot(x,Y3(:,i))
%     axis([0,1,0,1]);
    title(sprintf('DEIM MOR FEM Animation t= %0.3f',((i-1)*(Paras.t_end/Paras.t_n))))   
    
    
    F(i) = getframe;
    
   
end

% movie(F,2)

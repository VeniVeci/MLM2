function [Y,T,Time_Ode_solver]=Burger1D_FEM_DBC_MOR_DEIM_SolverF(Paras,U,U_DEIM,P)

% Burgers equation 1D case, finite element method, Dirichlet boundary
% conditions & homeogeneous B.C, solver.

% Problem model:
% [u(t,x)]_t+1/2[u^2(t,x)]_x-q [u(t,x)]_xx=f(t,x), x \in (0,1),t>0;

% Boundary Condition (Dirichlet):
% u(t,0)=0; u(t,1)=0;

% Initial Conditions:
% u(0,x)=u_0(x);

% close all
% clear;

%% Setup
% ------------------Problem Parameters------------------------------------- 
Re =Paras.Re;    % Reynolds Number
v=1/Re;     % viscosity

u0a=Paras.u0a; 
u0b=Paras.u0b; 
% ------------------Solver Parameters--------------------------------------
n=Paras.n;           % Total Spatial elements
t_end=Paras.t_end;        % End time
t_n=Paras.t_n;        % Number of time step 
t=0:(t_end/t_n):t_end; % time sequence (t=0 does not account 1 step)

% solver = 'ode45';
% options = odeset('RelTol',1e-6,'AbsTol',1e-10);

% ------------------Calculating temp variable----------------------------- 
h=1/n;      % space step size
x = 0:h:1;  % coordinate sequence

%% Main 

% Generate Mass matrix;
M=zeros(n-1,n-1);    
for i=1:n-1
    M(i,i)=2/3;
    if i>1
        M(i-1,i)=1/6;
    end
    if i<n-1
        M(i+1,i)=1/6;
    end
end
M=M*h;

% for i=1:n-1
%     for j=1:n-1
%         M(i,j)=integral(@(x)PhiPhi(x,i,j,n),0,1);
%     end
% end


% Generate B (Convective Term)
B=zeros(n-1,n-1);   
for i=1:n-1
    if i>1
        B(i-1,i)=1/2;
    end
    if i<n-1
        B(i+1,i)=-1/2;
    end
end

% for i=1:n-1
%     for j=1:n-1
%         B(i,j)=integral(@(x)PhiPhiDiff(x,i,j,n),0,1);
%     end
% end


% Generate Stiffness Matrix
C=zeros(n-1,n-1);
for i=1:n-1
    C(i,i)=2;
    if i>1
        C(i-1,i)=-1;
    end
    if i<n-1
        C(i+1,i)=-1;
    end
end
C=C*h^(-1);

% for i=1:n-1
%     for j=1:n-1
%         C(i,j)=integral(@(x)PhiDiffPhiDiff(x,i,j,n),0,1);
%     end
% end

% Generate the G vector
G=zeros(n-1,1);
for i=1:n-1
        G(i,1)=integral(@(x)PhiU0(x,i,n,u0a,u0b),0,1);
end

% Generate F matrix. The source term is assumed zero so F=0 for all time.
F=zeros(n-1,1);
for i=1:n-1
%     F(i,1)=integral(@(x)(0.02*exp(x).*FEM_HatFunc(i,n,x)),0,1);      
    F(i,1)=integral(@(x)(0.0*exp(x).*FEM_HatFunc(i,n,x)),0,1);  
end



% Generate initial condition
y_int=M\G;


% model reduction
Z_int=U'*y_int;

Mr=U'*M*U;
Fr=U'*F;
Cr=U'*C*U;
Br=U'*B;


%DEIM Dr
Dr=U_DEIM*inv(P'*U_DEIM);

% options = odeset('RelTol',1e-1,'AbsTol',1e-1);

tic;
[T,Z] = ode45(@(t,z) Burgers1D_DBC_FEM_ODE_MOR_DEIM_func(t,z,Mr,Br,Cr,Fr,v,U,Dr,P),t,Z_int); % Solve ODE
% [T,Z] = ode45(@(t,z) Burgers1D_DBC_FEM_ODE_MOR_DEIM_func(t,z,Mr,Br,Cr,Fr,v,U,Dr,P),t,Z_int,options); % Solve ODE
% [T,Z] = ode15s(@(t,z) Burgers1D_DBC_FEM_ODE_MOR_DEIM_func(t,z,Mr,Br,Cr,Fr,v,U,Dr,P),t,Z_int,options); % Solve ODE
% [T,Z] = ode23tb(@(t,z) Burgers1D_DBC_FEM_ODE_MOR_DEIM_func(t,z,Mr,Br,Cr,Fr,v,U,Dr,P),t,Z_int); % Solve ODE
Time_Ode_solver=toc;


Z=Z';
Y=U*Z;


% Add Dirichlet Boundary value
[nRow,ncolumn]=size(Y);
if (ncolumn~=(t_n+1)), Y=zeros(Paras.n-1,t_n+1);end


Y=[zeros(1,t_n+1);Y;zeros(1,t_n+1)];


% Ploting

% figure(1)
% meshc(Y(:,1:100))
% figure(1)
% meshc(Y)
% title(sprintf('FEM solution'))

% for i = 1:t_n+1
%     figure(3)
%     plot(x,Y(:,i))
%     axis([0,1,-3,3]);
%     title(sprintf('Animation t= %0.3f',((i-1)*(t_end/t_n))))
%     pause(0.1);  
% end



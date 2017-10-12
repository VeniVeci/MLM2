function [Y,T,Time_Ode_solver]=Burger1D_FEM_NDBC_SolverF(Paras)

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
% M=zeros(n-1,n-1);    
% for i=1:n-1
%     M(i,i)=2/3;
%     if i>1
%         M(i-1,i)=1/6;
%     end
%     if i<n-1
%         M(i+1,i)=1/6;
%     end
% end
% M=M*h;

for i=0:n-1
    for j=0:n-1
%         M(i,j)=integral(@(x)PhiPhi(x,i,j,n),0,1);
        M(i+1,j+1)=integral(@(x)(FEM_HatFunc(j,n,x).*FEM_HatFunc(i,n,x)),0,1);       
    end
end


% Generate B (Convective Term)
% B=zeros(n-1,n-1);   
% for i=1:n-1
%     if i>1
%         B(i-1,i)=1/2;
%     end
%     if i<n-1
%         B(i+1,i)=-1/2;
%     end
% end

for i=0:n-1
    for j=0:n-1
        B(i+1,j+1)=integral(@(x)(FEM_HatFunc(i,n,x ).*FEM_HatFunc_Diff(j,n,x )),0,1);
    end
end


% Generate Stiffness Matrix
% C=zeros(n-1,n-1);
% for i=1:n-1
%     C(i,i)=2;
%     if i>1
%         C(i-1,i)=-1;
%     end
%     if i<n-1
%         C(i+1,i)=-1;
%     end
% end
% C=C*h^(-1);

for i=0:n-1
    for j=0:n-1
        C(i+1,j+1)=integral(@(x)(FEM_HatFunc_Diff(i,n,x ).*FEM_HatFunc_Diff(j,n,x)),0,1);
    end
end

% Generate the G vector
G=zeros(n,1);
for i=0:n-1
%         G(i,1)=integral(@(x)PhiU0(x,i,n,u0a,u0b),0,1);
        G(i+1,1)=integral(@(x)(U0(x,u0a,u0b).*FEM_HatFunc(i,n,x)),0,1);       
end

% Generate F matrix. The source term is assumed zero so F=0 for all time.
% F=zeros(n-1,1);

for i=0:n-1
    F(i+1,1)=integral(@(x)(0.01*exp(x).*FEM_HatFunc(i,n,x)),0,1);        
end


D=zeros(n,1);
D(1,1)=-0.01;



% Generate initial condition
y_int=M\G;



tic;
[T,Y] = ode45(@(t,y) Burgers1D_NDBC_FEM_ODE_func(t,y,M,B,C,F,v,D),t,y_int); % Solve ODE
Time_Ode_solver=toc;

Y=Y';

% Add Dirichlet Boundary value
%tempy=[zeros(1,length(t));y;zeros(1,length(t))];


Y=[Y;zeros(1,t_n+1)];



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



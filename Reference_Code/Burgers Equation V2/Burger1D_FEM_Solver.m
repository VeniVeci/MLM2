function [Y,T,Time]=Burger1D_FEM_Solver(Paras)

% Burgers equation 1D case, finite element method, Dirichlet & Neumann 
% boundary conditions solver.

% Problem model:
% [u(t,x)]_t+1/2[u^2(t,x)]_x-q [u(t,x)]_xx=f(t,x), x \in (0,1),t>0;

% Boundary Condition (Dirichlet):
% u(t,0)=Paras.BC.ux0; u(t,1)=Paras.BC.ux1;
% u(t,0)=Paras.BC.dx_ux0; u(t,1)=Paras.BC.dx_ux1;

% Initial Conditions:
% u(0,x)=Paras.u_0(x);

% close all

%% Setup
% ------------------Problem Parameters------------------------------------- 
Re =Paras.Re;    % Reynolds Number
v=1/Re;     % viscosity

u0a=Paras.u0a; 
u0b=Paras.u0b; 

Sa=Paras.Sa;
Sb=Paras.Sb;

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

M=zeros(n+1,n+1); 
for i=0:n
    for j=0:n
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

B=zeros(n+1,n+1);   
for i=0:n
    for j=0:n
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

C=zeros(n+1,n+1);
for i=0:n
    for j=0:n
        C(i+1,j+1)=integral(@(x)(FEM_HatFunc_Diff(i,n,x ).*FEM_HatFunc_Diff(j,n,x)),0,1);
    end
end

% Generate the G vector
g=zeros(n+1,1);
for i=0:n
%         G(i,1)=integral(@(x)PhiU0(x,i,n,u0a,u0b),0,1);
        g(i+1,1)=integral(@(x)(U0(x,u0a,u0b).*FEM_HatFunc(i,n,x)),0,1);       
end

% Generate sources F matrix. The source term is assumed zero so F=0 for all time.

 F=zeros(n+1,1);
for i=0:n-1
    F(i+1,1)=integral(@(x)(Sources(x,Paras.Sa,Paras.Sb).*FEM_HatFunc(i,n,x)),0,1);        
end


% Nuemann source term d vector
d=zeros(n+1,1);
d(1,1)=  Paras.BC.dx_ux0;
d(n+1,1)=Paras.BC.dx_ux1;

if Paras.BC.dx_ux0==0, nStart=2; else nStart=1; end; %If Dirichlet BC   
if Paras.BC.dx_ux1==0, nEnd=n; else nEnd=n+1;   end; %If Dirichlet BC

y_int=M(nStart:nEnd,nStart:nEnd)\g(nStart:nEnd,1);
% y_int2=M\g;
% U=zeros(n+1,t_n);

tic;

% [T,Y] = ode45(@(t,y) Burgers1D_NDBC_FEM_ODE_func(t,y,M(nStart:nEnd,nStart:nEnd),B(nStart:nEnd,nStart:nEnd),C(nStart:nEnd,nStart:nEnd),F(nStart:nEnd,1),v,d(nStart:nEnd,1)),t,y_int); % Solve ODE
[T,Y] = ode45(@(t,y) Burgers1D_ODE_func(t,y,M(nStart:nEnd,nStart:nEnd),B(nStart:nEnd,nStart:nEnd),C(nStart:nEnd,nStart:nEnd),F(nStart:nEnd,1),v,d(nStart:nEnd,1)),t,y_int); % Solve ODE
Time=toc;
Y=Y';

if Paras.BC.dx_ux0==0
    Y=[Paras.BC.ux0.*ones(1,t_n+1);Y];
end
if Paras.BC.dx_ux1==0
    Y=[Y;Paras.BC.ux1.*ones(1,t_n+1)];
end


end

function [dy]=Burgers1D_ODE_func(t,y,M,B,C,F,v,D)
    dy = M\ (F-(1/2)*B*y.^2-v*(D+C*y));
end

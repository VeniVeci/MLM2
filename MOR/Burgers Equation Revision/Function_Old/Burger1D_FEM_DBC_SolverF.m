function [Y,T,Time_Ode_solver]=Burger1D_FEM_DBC_SolverF(Paras,gx,u0)
%%
% Burgers equation 1D case, finite element method, Dirichlet boundary
% conditions & homeogeneous B.C, solver.
% 
% Problem model:
% [u(t,x)]_t+1/2[u^2(t,x)]_x-q [u(t,x)]_xx=f(t,x), x \in (0,1),t>0;
% 
% Boundary Condition (Dirichlet):
% u(t,0)=0; u(t,1)=0;
% 
% Initial Conditions:
% u(0,x)=u_0(x);
% 
% 
%
% History:
% Lost date,  WeiX, first edition 
% 10-01-2017, WeiX, modify structure


%% Setup
% ------------------Problem Parameters------------------------------------- 
Re =Paras.Re;    % Reynolds Number
v=1/Re;     % viscosity

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
%         G(i,1)=integral(@(x)PhiU0(x,i,n,u0a,u0b),0,1);
        G(i,1)=integral(@(x)(u0(x).*FEM_HatFunc(i,n,x)),0,1); 
end

% Generate F matrix. The source term is g(x) a function
F=zeros(n-1,1);
for i=1:n-1
%     F(i,1)=integral(@(x)(0.01*exp(x).*FEM_HatFunc(i,n,x)),0,1);    
    F(i,1)=integral(@(x)(gx(x).*FEM_HatFunc(i,n,x)),0,1);
end

% Generate initial condition
y_int=M\G;



tic;
[T,Y] = ode45(@(t,y) Burgers1D_DBC_FEM_ODE_func(t,y,M,B,C,F,v),t,y_int); % Solve ODE
Time_Ode_solver=toc;

Y=Y';

% Add Dirichlet Boundary value
%tempy=[zeros(1,length(t));y;zeros(1,length(t))];
Y=[zeros(1,t_n+1);Y;zeros(1,t_n+1)];

end 




function  Phi  = FEM_HatFunc(i,N,x )
% Description:
% Finite element widely used hat function as the basic function3
% a better way to express FEM hat function.(simpler, without 'if' operation
% and other logical operation, Thus can be used with ingegral function)
%
% WARMING: inputs should always satisfy: x>=0 && x<=1; 
%
% Synopsis:
% Phi  = FEM_HatFunc3(i,N,x )
%
% Input: 
%       i                 % index of the ith basic subfunction  i \in [0,N];
%       N                 % number of element. Total number of basic 
%                         % function =N+1.(including i=0) 
%                         %(normally N=number of element+1 in 1D case) 
%       x                 % variable x
% Output: 
%       Phi               % output 
% Pcakage Require:
%
% Modifications:
% 27-May-2015, WeiX, first edition 

h=1/N;

Phi=1-abs(x-i*h)/h;

% for matrix computation
Phi=Phi.*(x<h*(i+1));
Phi=Phi.*(x>h*(i-1));

%Phi=Phi.*(heaviside(h*(i-1))-heaviside(h*(i+1))) % A alternative solution.
%But wrong. Need improve.

% Equivalent to the following but capible for matrix computation.
% h=1/N;
% X=h:h:1;
% 
% if i ==0
%   
%     if ((x>=0) && (x<=X(1)))
%         Phi=(X(1)-x)/h;
%     else 
%         Phi=0;
%     end
%     
% elseif i==1    %Compromise for X(0)(=0) not access in MATLAB
%     if ((x>=0) && (x<=X(i)))   
%         Phi=(x-0)/h;    
%     elseif x>=X(i) && x<=X(i+1)
%         Phi=(X(i+1)-x)/h;
%     else
%         Phi=0;        
%     end        
%        
% elseif i>=1 && i<=N-1
%     
%     if x>X(i-1) && x<=X(i)   
%         Phi=(x-X(i-1))/h;     
%         
%     elseif x>=X(i) && x<=X(i+1)
%         Phi=(X(i+1)-x)/h;
%     else
%         Phi=0;        
%     end    
%     
% elseif i ==N
%     
%     if x>=X(N-1) && x<=X(N)
%         Phi=(x-X(N-1))/h;
%     else 
%         Phi=0;
%     end
% 
% end

end


function [dy]=Burgers1D_DBC_FEM_ODE_func(t,y,M,B,C,F,v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Burgers' equation 1D, finite element method, Dirichlet boundary 
% conditions, solver.
% 
% Function inputs:
%
% t : current time
% y : current y
% v : viscosity
%
%
% Function outputs:
% out y : computed derivative
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %F = RHS F(f, in t);
    % For the moment assume F=0;

    dy = M\ (F-(1/2)*B*y.^2-v*C*y);

end


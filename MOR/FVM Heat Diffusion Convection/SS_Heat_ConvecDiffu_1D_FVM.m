%Simulating the steady state heat Convection-Diffusion 1D equation  
...by the Finite Volume Method
    
% Governning equation : k*Txx+S=(rho*u*phi)x; Txx=dT/d^2(x);
...d(rho*u*phi)/dx
    
% The simulation would be stable when F/D=Pe<2.(Peclet number)
...Otherwise it will ocilate. An upwind strategy is required.

clear

% -----------------------Setup Parameters-------------------------------
L=1;        % total length is 0.1m
Area=0.01;  % Cross-session area is 0.01m^2.
k=0.1;      % Thermal conductivity is 1000W/m/K
T_L=1;    % Left boundary temperature
T_R=0;    % Right boundary temperature
u=0.1;      % Constant flow velosity (from left to right)
rho=1;      % flow density

Num_node=50;      % number of volumes

% -------Initializing Parameters------
dx=L/Num_node;
F=rho*u;
D=k/dx;

% -----------------------Initializing -------------------------------
X=zeros(Num_node,1);
A=zeros(Num_node,Num_node);
b=zeros(Num_node,1);

% -----------------------F V M------------------------------------------
% Problem to be solved is vector X. It is solved by AX=b. A,b are known.
%--------------------Build Matrix A-------------------------------------
%---For iner node-------   
for i = 2:Num_node-1
    
    a_W=D+F/2;
    a_E=D-F/2;
    a_P=-a_W-a_E;
    
    A(i,i-1)=a_W;
    A(i,i)  =a_P;
    A(i,i+1)=a_E;
end
%---For Boundary node-----------------------------
%---For Left node---------------------
    a_W=0;
    a_E=D-F/2;
    S_P=-(2*D+F);
    a_P=S_P-a_W-a_E;
    %S_u=2*k*A/dx*T_A
    
    i=1;
   %A(i,i-1)=a_W;
    A(i,i)  =a_P;
    A(i,i+1)=a_E;
%---For Right node---------------------
    a_W=D+F/2;
    a_E=0;
    S_P=-(2*D-F);
    a_P=S_P-a_W-a_E;
    %S_u=2*k*A/dx*T_B
    
    i=Num_node;
    A(i,i-1)=a_W;
    A(i,i)  =a_P;
   %A(i,i+1)=a_E;
   
%--------------------Build Vector b-------------------------------------
%---For Left node---------------------
S_u=(2*D+F)*T_L;
i=1;
b(i,1)=-S_u;
%---For Right node---------------------
S_u=(2*D-F)*T_R;
i=Num_node;
b(i,1)=-S_u;

%-------------------Solve AX=b----------------------------------------
% X=inv(A)*b;
X=A\b;
%-------------------Plot----------------------------------------------
X_label=(dx/2:dx:L-dx/2)';
plot(X_label,X,'-o');

%------------------Exact Solution------------------------------------
x=X_label;
a=exp(rho*u*x/k)-1;
b=exp(rho*u*L/k)-1
T=(a./b)*(T_R-T_L)+T_L;
hold on 
plot(X_label,T,'-r+');
hold off
%------------------End-------------------------------------------------







    
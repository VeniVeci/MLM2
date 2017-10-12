%Simulating the steady state 1-D Diffusion equation (Fourier's equation)
...by the Finite Volume Method
    
% Governning equation : k*Txx+S=0; Txx=dT/d^2(x)

clear
%%
% -----------------------Setup Parameters-------------------------------
L=1;      % total length is 1m
Area=0.01;  % Cross-session area is 0.01m^2.
k=1000;     % Thermal conductivity is 1000W/m/K
T_L=100;    % Left boundary temperature
T_R=500;    % Right boundary temperature

Num_node=5;      % number of volumes

% -----------------------Initializing -------------------------------
X=zeros(Num_node,1);
A=zeros(Num_node,Num_node);
b=zeros(Num_node,1);
dx=L/Num_node;

%%
% -----------------------F V M------------------------------------------
% Problem to be solved is vector X. It is solved by AX=b. A,b are known.
%--------------------Build Matrix A-------------------------------------
%---For iner node-------   
for i = 2:Num_node-1
    
    a_W=k*Area/dx;
    a_E=k*Area/dx;
    a_P=-a_W-a_E;
    
    A(i,i-1)=a_W;
    A(i,i)  =a_P;
    A(i,i+1)=a_E;
end
%---For Boundary node-----------------------------
%---For Left node---------------------
    a_W=0;
    a_E=k*Area/dx;
    S_P=-2*k*Area/dx;
    a_P=S_P-a_W-a_E;
    %S_u=2*k*A/dx*T_A
    
    i=1;
   %A(i,i-1)=a_W;
    A(i,i)  =a_P;
    A(i,i+1)=a_E;
%---For Right node---------------------
    a_W=k*Area/dx;
    a_E=0;
    S_P=-2*k*Area/dx;
    a_P=S_P-a_W-a_E;
    %S_u=2*k*A/dx*T_B
    
    i=Num_node;
    A(i,i-1)=a_W;
    A(i,i)  =a_P;
   %A(i,i+1)=a_E;
   
%--------------------Build Vector b-------------------------------------
%---For Left node---------------------
S_u=2*k*Area/dx*T_L;
i=1;
b(i,1)=-S_u;
%---For Right node---------------------
S_u=2*k*Area/dx*T_R;
i=Num_node;
b(i,1)=-S_u;

%-------------------Solve AX=b----------------------------------------
X=inv(A)*b;

%%
%-------------------Plot----------------------------------------------
X_label=(dx/2:dx:L-dx/2)';
plot(X_label,X,'-o');

    
%Simulating the unsteady state 1-D Diffusion equation (Fourier's equation)
...by the Finite Volume Method using fully implicit scheme
    
% Governning equation : k*Txx+S=0; Txx=dT/d^2(x)

clear
%%
% ------Setup problem Parameters-------------------------------
L=1;      % total length is 1m
Area=0.01;  % Cross-session area is 0.01m^2.
k=1000;     % Thermal conductivity is 1000W/m/K
rhoc=1e6;     % density times conductivity
T_L=100;    % Left boundary temperature
T_R=500;    % Right boundary temperature

t_end=300;       % total simulation time 
dt=2;   % time step

T_0=1;      % initial temperature

% ------Setup solver Parameters---------------------------------
Num_node=1000;      % number of volumes

% -----------------------Initializing -------------------------------
X=zeros(Num_node,1);
A=zeros(Num_node,Num_node);
b=zeros(Num_node,1);
dx=L/Num_node;

X_0=ones(Num_node,1)*T_0;

%%

Index_frame=1;
for t=0:dt:t_end

    % -----------------------F V M------------------------------------------
    % Problem to be solved is vector X. It is solved by AX=b. A,b are known.
    %--------------------Build Matrix A-------------------------------------
    %---For iner node-------   
    for i = 2:Num_node-1

        a_W=k/dx;
        a_E=k/dx;
        
        a_P0=rhoc*dx/dt;
        a_P=-a_P0-a_W-a_E;

        A(i,i-1)=a_W;
        A(i,i)  =a_P;
        A(i,i+1)=a_E;
    end
    %---For Boundary node-----------------------------
    %---For Left node---------------------
        a_W=0;
        a_E=k/dx;
        S_P=-2*k/dx;
        a_P0=rhoc*dx/dt;
        a_P=S_P-a_P0-a_W-a_E;
        %S_u=2*k*A/dx*T_A

        i=1;
       %A(i,i-1)=a_W;
        A(i,i)  =a_P;
        A(i,i+1)=a_E;
    %---For Right node---------------------
        a_W=k/dx;
        a_E=0;
        S_P=-2*k/dx;
        a_P0=rhoc*dx/dt;
        a_P=S_P-a_P0-a_W-a_E;
        %S_u=2*k*A/dx*T_B

        i=Num_node;
        A(i,i-1)=a_W;
        A(i,i)  =a_P;
        %A(i,i+1)=a_E;

    %--------------------Build Vector b-------------------------------------
    a_P0=rhoc*dx/dt;
    for i = 2:Num_node-1
        b(i,1)=-a_P0*X_0(i);
    end
    %---For Left node---------------------
    S_u=2*k/dx*T_L;
    i=1;
    b(i,1)=-S_u-a_P0*X_0(i);
    %---For Right node---------------------
    S_u=2*k/dx*T_R;
    i=Num_node;
    b(i,1)=-S_u-a_P0*X_0(i);
       
    %-------------------Solve AX=b----------------------------------------
    X=A\b;
    
    X_0=X;
    %%
    %-------------------Plot----------------------------------------------
    X_label=(dx/2:dx:L-dx/2)';
    h=plot(X_label,X,'-o');
%     pause(0.1)
     F(Index_frame)=getframe;
     Index_frame=Index_frame+1;
    

end
    pause(1)
    movie(F)
    
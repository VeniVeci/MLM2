%Simulating the unsteady state 1-D Diffusion equation (Fourier's equation) 
...by the Finite Volume Method using fully implicit scheme. (1200 times faster for pure calculation( without ploting)!!!)
    
% Governning equation : k*Txx+S=0; Txx=dT/d^2(x)

clear
%%
% ------Setup problem Parameters-------------------------------
L=1;        % total length is 1m
Area=0.01;  % Cross-session area is 0.01m^2.
k=1000;     % Thermal conductivity is 1000W/m/K
rhoc=1e6;   % density times conductivity
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

X_int=ones(Num_node,1)*T_0;

%%

start_time = cputime;
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
  
    SU=zeros(Num_node,1);
    %---For Left node---------------------
    S_U(1,1)=2*k/dx*T_L;
    %---For Right node---------------------
    S_U(Num_node,1)=2*k/dx*T_R;
       
    
Index_frame=1;
X_0=X_int;
Rec_X=X_0;
for t=0:dt:t_end
    
    b=-S_U-a_P0*X_0;
    %-------------------Solve AX=b----------------------------------------
    X=A\b;
    
    Rec_X=[Rec_X,X];
    X_0=X;
    %%
%     %-------------------Plot----------------------------------------------
%     X_label=(dx/2:dx:L-dx/2)';
%     h=plot(X_label,X,'-o');
% %     pause(0.1)
%      F(Index_frame)=getframe;
%      Index_frame=Index_frame+1;
    
end
Time_Normal=cputime-start_time




%% MOR mode 
approximate_degree=10;

% -----------------------Initializing -------------------------------
X_0=ones(Num_node,1)*T_0;


[U,S,~]=svd(Rec_X);  % U*S*V'=Rec_X
U=U(:,1:approximate_degree);
eigenvalues=diag(S);



start_time = cputime;
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
  
    SU=zeros(Num_node,1);
    %---For Left node---------------------
    S_U(1,1)=2*k/dx*T_L;
    %---For Right node---------------------
    S_U(Num_node,1)=2*k/dx*T_R;
             
    %------------------MOR---------------------------------
    A_Mor=U'*A*U;
    S_U_Mor=U'*S_U;
    
%     X_Mor_0=U'*X_0;
    
    
Index_frame=1;
X_0=X_int;
X_Mor_0=U'*X_0;
Rec_X_Mor=X_Mor_0;
for t=0:dt:t_end    
    
    %-------------------Solve AX=b----------------------------------------
    
%     b=-S_U-a_P0*X_0;
%     X=A\b;
%     Rec_X=[Rec_X,X_0];
%     X_0=X;
    
%     X_Mor_0=U'*X_Mor;
%     b_Mor=-U'*S_U-a_P0*X_Mor_0;
    
    
    b_Mor=-S_U_Mor-a_P0*X_Mor_0;
    X_Mor=A_Mor\b_Mor;
    Rec_X_Mor=[Rec_X_Mor,X_Mor];
    X_Mor_0=X_Mor;
    
    
%     b_Mor=U'*b;
%     X_Mor=A_Mor\b_Mor;
%     
%     X=U*X_Mor;
%     
%     Rec_X_Mor=[Rec_X_Mor,X_0];
%     X_0=X;
    %%
    
%     %-------------------Plot----------------------------------------------
%     X=U*X_Mor;
%     X_label=(dx/2:dx:L-dx/2)';
%     h=plot(X_label,X,'-o');
% %     pause(0.1)
%      F_Mor(Index_frame)=getframe;
%      Index_frame=Index_frame+1;
    
end
Rec_X_byMor=U*Rec_X_Mor;

Time_MOR=cputime-start_time


%     pause(1)
%     movie(F)
    
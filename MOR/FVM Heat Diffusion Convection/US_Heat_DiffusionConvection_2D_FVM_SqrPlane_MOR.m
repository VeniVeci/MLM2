%Simulating the unsteady state heat flow Diffusion Convection 2D equation  
...by the Finite Volume Method using centreal differencing scheme in space
...and fully implicit scheme in time.

% Instruction: Remove/comment the plot code would should the Huge different
% between ROM and HDM


    
% Governning equation : d(rho*T)/dt+k*Txx+k*Tyy=(rho*u*phi)x+(rho*u*phi)y+S; 
...Txx=dT/d^2(x);Tyy=dT/d^2(y);d(rho*u*phi)/dx;d(rho*u*phi)/dy;
...u=(u,y)' a vector indicates velocity.

% The point are discritized as follow grid. e.g a 3*4 grid. Numbers are
% index of point
% 
%        4--8--12
%        |  |   | 
%        3--7--11
%        |  |  |
%     ^  2--6--10
%     |  |  |  |
%  u_y|  1--5--9
%        ----> u_x
%
% Pe_x=F_x/D_x=(rho*u_x*dy)/(k*dy/dx)=rho*u_x*dx/k 
... Only if Pe_x<2 the result is stable and accurate in x direction
% Pe_y=F_y/D_y=(rho*u_y*dx)/(k*dx/dy)=rho*u_y*dy/k
... Only if Pe_y<2 the result is stable and accurate in y direction
%     
% Modifications:
% 1-April-2014, WeiX, first edition 

% close all
clear

%% -----------------------Setup Parameters----------------------------------
%  -----------------------Problem Parameters---------------------------
Lx=10;                % length of x axis is 1m
Ly=10;                % length of y axis is 1m
k=1;               % Thermal conductivity is 1000W/m/K
rho=1;                 % flow density
T_L=100;               % Temperature of Left    boundary 
T_R=400;               % Temperature of Right   boundary 
T_B=100;               % Temperature of Bottom  boundary 
T_T=400;               % Temperature of Top     boundary 
u_x=1;                 % Velocity x direction 
u_y=-1;                 % Velocity y direction 

T_0=1000;              % initial temperature

t_end=15;             % total simulation time 
dt=0.1;                  % time step

% -----------------------Solver Parameters---------------------------
Num_node_x=20;        % number of volumes on x direction
Num_node_y=20;        % number of volumes on y direction

%% -------Initializing computational Parameters------
dx=Lx/Num_node_x;
dy=Ly/Num_node_y;
Num_node_total=Num_node_x*Num_node_y;

% --------Initializing memory -------------------------------
X=zeros(Num_node_total,1);
X_int=rand(Num_node_total,1)*T_0;
A=zeros(Num_node_total,Num_node_total);
b=zeros(Num_node_total,1);

% --------PreProcessing -------------------------------------
Pe_x=(rho*u_x*dy)/(k*dy/dx);
Pe_y=(rho*u_y*dx)/(k*dx/dy);
if Pe_x>2 warning('Pe_x=%d >2, result of x-direction is untable and non accurate',Pe_x); end
if Pe_y>2 warning('Pe_y=%d >2, result of y-direction is untable non non accurate',Pe_y); end
    
%% 
Index_frame=1;

    % -----------------------F V M Solver-------------------------------------
    % Problem to be solved is vector X. It is solved by AX=b. A,b are known.
    %
    %--------------------Build Matrix A-------------------------------------
    Node_P_index=0;
    for i = 1:Num_node_x
        for j = 1:Num_node_y

            Node_P_index=Node_P_index+1;            % point index 

            P(1,Node_P_index)=dx*(1/2+(i-1));       % P(1,:) is the X cooridnate of points
            P(2,Node_P_index)=dy*(1/2+(j-1));       % P(2,:) is the y cooridnate of points
            P(3,Node_P_index)=Node_P_index;         % P(3,:) is the index of points

            Node_W_index=Node_P_index-Num_node_y;   % identify west  neighbour (in a retangular grid)
            Node_E_index=Node_P_index+Num_node_y;   % identify east  neighbour (in a retangular grid)
            Node_S_index=Node_P_index-1;            % identify south neighbour (in a retangular grid)
            Node_N_index=Node_P_index+1;            % identify north neighbour (in a retangular grid)

            F_w=rho*u_x*dy;
            F_e=rho*u_x*dy;
            F_s=rho*u_y*dx;
            F_n=rho*u_y*dx;

            D_w=k*dy/dx;                            
            D_e=k*dy/dx;
            D_s=k*dx/dy;
            D_n=k*dx/dy;

            a_W=D_w+F_w/2;                          % set default coefficients which is used for a inner node.
            a_E=D_e-F_w/2;
            a_S=D_s+F_s/2;
            a_N=D_n-F_n/2;

            S_PW=0;
            S_PE=0;
            S_PS=0;
            S_PN=0;

            if i==1                     % point on the left (next to boundary)
                a_W=0;                  % set coefficient a_W zero
                S_PW=-2*D_w-F_w;        % sources form the west boundary 
                Node_W_index=0;         % west neighbour is void 
            end          

            if i==Num_node_x            % point on the right
                a_E=0;
                S_PE=-2*D_e+F_e; 
                Node_E_index=0;      
            end         

            if j==1                     % point on the bottom
                a_S=0;
                S_PS=-2*D_s-F_s;
                Node_S_index=0;      
            end         

            if j==Num_node_y            % point on the top
                a_N=0;
                S_PN=-2*D_n+F_n;
                Node_N_index=0;      
            end          
            
            a_P0=rho*dx*dy/dt;
            
            S_P=S_PW+S_PE+S_PS+S_PN;    
            S_u=-S_PW*T_L-S_PE*T_R-S_PS*T_B-S_PN*T_T;
            a_P=S_P-a_P0-a_W-a_E-a_S-a_N+F_w-F_e+F_s-F_n;

           A(Node_P_index,Node_P_index)=a_P;        % putting coefficients in matrix A
           if Node_W_index~=0   A(Node_P_index,Node_W_index)=a_W;end
           if Node_E_index~=0   A(Node_P_index,Node_E_index)=a_E;end
           if Node_S_index~=0   A(Node_P_index,Node_S_index)=a_S;end
           if Node_N_index~=0   A(Node_P_index,Node_N_index)=a_N;end

%            b(Node_P_index,1)=-S_U-a_P0*X_0(Node_P_index);                  % putting coefficients in matrix b
           
           S_U(Node_P_index,1)=S_u;

        end
    end
    
    
X_0=X_int;
Rec_X=X_0;

Index_frame=1;
start_time = cputime;
    
for t=0:dt:t_end
    
    b=-S_U-a_P0*X_0; 
    %-------------------Solve AX=b----------------------------------------
    % X=inv(A)*b;
    X=A\b;
    % X=linsolve(A,b);    % For better performance
    Rec_X=[Rec_X,X];
    X_0=X;

    %%
    %-------------------Plot----------------------------------------------

    Tx=reshape(P(1,:),Num_node_y,Num_node_x);       % Rearrange coordinate
    Ty=reshape(P(2,:),Num_node_y,Num_node_x);
    T =reshape(X,Num_node_y,Num_node_x);

    figure(1)
    subplot(2,1,1),contour(Tx,Ty,T),
    title('Temperature (Steady State)'),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
    subplot(2,1,2),pcolor(Tx,Ty,T),shading interp,
    title('Temperature (Steady State)'),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal

%     figure(2)
%     pcolor(Tx,Ty,T),shading interp,
%     title('Temperature (Steady State)'),xlabel('x'),ylabel('y'),colorbar,axis([0,Lx,0,Ly]),axis equal;
% 
%     % figure(3)
%     % mesh(Tx,Ty,T),shading interp,
%     figure(4)
%     surf(Tx,Ty,T),shading interp,
%     figure(5)
%     contour(Tx,Ty,T)


     F(Index_frame)=getframe;
     Index_frame=Index_frame+1;
    t
    
end
Time_Normal=cputime-start_time
%     pause(1)
%     movie(F)




%% MOR 
approximate_degree=10;

% [~,series]=size(Rec_X);
% Num_snapshot=30;


[U,S,~]=svd(Rec_X);  % U*S*V'=Rec_X
U=U(:,1:approximate_degree);
eigenvalues=diag(S);


%------------------MOR---------------------------------
A_Mor=U'*A*U;
S_U_Mor=U'*S_U;


X_0=X_int;
X_Mor_0=U'*X_0;
Rec_X_Mor=X_Mor_0;

Index_frame=1;
start_time = cputime;
for t=0:dt:t_end
   
    
    b_Mor=-S_U_Mor-a_P0*X_Mor_0;
    X_Mor=A_Mor\b_Mor;
    Rec_X_Mor=[Rec_X_Mor,X_Mor];
    X_Mor_0=X_Mor;
    
    %%
    %-------------------Plot----------------------------------------------
    X=U*X_Mor;
    
    Tx=reshape(P(1,:),Num_node_y,Num_node_x);       % Rearrange coordinate
    Ty=reshape(P(2,:),Num_node_y,Num_node_x);
    T =reshape(X,Num_node_y,Num_node_x);

    figure(1)
    subplot(2,1,1),contour(Tx,Ty,T),
    title('Temperature (Steady State)'),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
    subplot(2,1,2),pcolor(Tx,Ty,T),shading interp,
    title('Temperature (Steady State)'),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal

%     figure(2)
%     pcolor(Tx,Ty,T),shading interp,
%     title('Temperature (Steady State)'),xlabel('x'),ylabel('y'),colorbar,axis([0,Lx,0,Ly]),axis equal;
% 
%     % figure(3)
%     % mesh(Tx,Ty,T),shading interp,
%     figure(4)
%     surf(Tx,Ty,T),shading interp,
%     figure(5)
%     contour(Tx,Ty,T)


     F(Index_frame)=getframe;
     Index_frame=Index_frame+1;
    t
   
    
    
end
Rec_X_byMor=U*Rec_X_Mor;

Time_MOR=cputime-start_time

SSE_X=(Rec_X_byMor-Rec_X).^2;

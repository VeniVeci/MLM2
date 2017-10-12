%Simulating the steady state heat Diffusion 2D equation  
...by the Finite Volume Method
    
% Governning equation : k*Txx+k*Tyy+S=0; Txx=dT/d^2(x);Tyy=dT/d^2(y);

% The point are discritized as follow grid. e.g a 3*4 grid. Numbers are
% index of point
% 
%     4--8--12
%     |  |   | 
%     3--7--11
%     |  |  |
%     2--6--10
%     |  |  |
%     1--5--9

clear
%%
% -----------------------Setup Parameters----------------------------------
% -----------------------Problem Parameters---------------------------
Lx=1;                % length of x axis is 1m
Ly=1;                % length of y axis is 1m
k=1000;               % Thermal conductivity is 1000W/m/K
T_L=200;               % Temperature of Left    boundary 
T_R=400;               % Temperature of Right   boundary 
T_B=100;               % Temperature of Bottom  boundary 
T_T=100;               % Temperature of Top     boundary 

% -----------------------Solver Parameters---------------------------
Num_node_x=10;        % number of volumes on x direction
Num_node_y=10;        % number of volumes on y direction

%%
% -------Initializing Parameters------
dx=Lx/Num_node_x;
dy=Ly/Num_node_y;
Num_node_total=Num_node_x*Num_node_y;

% -----------------------Initializing -------------------------------
X=zeros(Num_node_total,1);
A=zeros(Num_node_total,Num_node_total);
b=zeros(Num_node_total,1);

% -----------------------F V M------------------------------------------
% Problem to be solved is vector X. It is solved by AX=b. A,b are known.
%--------------------Build Matrix A-------------------------------------
Node_P_index=0;
for i = 1:Num_node_x
    for j = 1:Num_node_y
                
        Node_P_index=Node_P_index+1;        % point index 
        
        P(1,Node_P_index)=dx*(1/2+(i-1));   % P(1,:) is the X cooridnate of points
        P(2,Node_P_index)=dy*(1/2+(j-1));   % P(2,:) is the y cooridnate of points
        P(3,Node_P_index)=Node_P_index;     % P(3,:) is the index of points
        
        Node_W_index=Node_P_index-Num_node_y;   % identify west  neighbour (in a retangular grid)
        Node_E_index=Node_P_index+Num_node_y;   % identify east  neighbour (in a retangular grid)
        Node_S_index=Node_P_index-1;            % identify south neighbour (in a retangular grid)
        Node_N_index=Node_P_index+1;            % identify north neighbour (in a retangular grid)

        a_W=k*dy/dx;                        % set default coefficients which is used for a inner node.
        a_E=k*dy/dx;
        a_S=k*dx/dy;
        a_N=k*dx/dy;
        
        S_PW=0;
        S_PE=0;
        S_PS=0;
        S_PN=0;
        
        if i==1                     % point on the left (next to boundary)
            a_W=0;                  % set coefficient a_W zero
            S_PW=-2*k*dy/dx;        % sources form the west boundary
            Node_W_index=0;         % west neighbour is void 
        end          
        
        if i==Num_node_x            % point on the right
            a_E=0;
            S_PE=-2*k*dy/dx; 
            Node_E_index=0;      
        end         
        
        if j==1                     % point on the bottom
            a_S=0;
            S_PS=-2*k*dx/dy;
            Node_S_index=0;      
        end         
        
        if j==Num_node_y            % point on the top
            a_N=0;
            S_PN=-2*k*dx/dy;
            Node_N_index=0;      
        end          
        
        S_P=S_PW+S_PE+S_PS+S_PN;    
        S_U=-S_PW*T_L-S_PE*T_R-S_PS*T_B-S_PN*T_T;
        a_P=S_P-a_W-a_E-a_S-a_N;

       A(Node_P_index,Node_P_index)=a_P;        % putting coefficients in matrix A
       if Node_W_index~=0   A(Node_P_index,Node_W_index)=a_W;end
       if Node_E_index~=0   A(Node_P_index,Node_E_index)=a_E;end
       if Node_S_index~=0   A(Node_P_index,Node_S_index)=a_S;end
       if Node_N_index~=0   A(Node_P_index,Node_N_index)=a_N;end
       
       b(Node_P_index,1)=-S_U;                  % putting coefficients in matrix b
                
    end
end

%-------------------Solve AX=b----------------------------------------
% X=inv(A)*b;
X=A\b;
% X=linsolve(A,b);    % For better performance


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

figure(2)
pcolor(Tx,Ty,T),shading interp,
title('Temperature (Steady State)'),xlabel('x'),ylabel('y'),colorbar,axis([0,Lx,0,Ly]),axis equal;

%%
%------------------Exact Solution------------------------------------

%------------------End-------------------------------------------------



%%
%------------------Backup code-----------------------------------------
% Node_P_index=0;
% for i = 1:Num_node_x
%     for j = 1:Num_node_y
%         
%        a_W=k*A/dx;
%        a_E=k*A/dx;
%        a_S=k*A/dy;
%        a_N=k*A/dy;
%         
%         if i==1                  % point on the left
%             if j==1              % point on the left bottom corner
%                 a_W=0;
%                 a_S=0;                          
%             elseif j==Num_node_y % point on the left top corner
%                 
%             else
%             end
%             
%         elseif i== Num_node_x    % point on the right
%             if j==1              % point on the right bottom corner
%             elseif j==Num_node_y % point on the top top corner
%             else
%             end
%             
%         elseif j==1              % point on the bottom
%             if i==1              % point on the bottom left corner
%             elseif j==Num_node_y % point on the top right corner
%             else
%             end
%             
%         elseif j== Num_node_y    % point on the top
%             if j==1              % point on the top left corner
%             elseif j==Num_node_y % point on the top right corner
%             else
%             end      
%         end     
%         Node_index=Node_index+1;
%     end
% end



    
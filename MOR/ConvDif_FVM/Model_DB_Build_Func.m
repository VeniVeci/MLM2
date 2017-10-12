function [ model ] = Model_DB_Build_Func( model,Paras,uParas)
% Dirichlay Boundary Condition
% Using grid infomation and model parameter to build model 
%
% Modifications:
% 6-April-2016, WeiX, first edition 

%% NOT FINISH.  DRAFT VERSION 

Index_frame=1;


for i=1:Paras.Num_node_total
    
    model.
    
    [u_wx,u_wy] = uxy(model-dx/2,P(2,Node_P_index),uParas);
    
    
    
    [u_wx,u_wy] = uxy(P(1,Node_P_index)-dx/2,P(2,Node_P_index),uParas);
    [u_ex,u_ey] = uxy(P(1,Node_P_index)+dx/2,P(2,Node_P_index),uParas);
    [u_sx,u_sy] = uxy(P(1,Node_P_index),P(2,Node_P_index)-dy/2,uParas);
    [u_nx,u_ny] = uxy(P(1,Node_P_index),P(2,Node_P_index)+dy/2,uParas);
    
        
   
    
    
    

    
    
end
    
    
    
    
    
    
    
    
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
            
            
            [u_wx,u_wy] = uxy(P(1,Node_P_index)-dx/2,P(2,Node_P_index),uParas);
            [u_ex,u_ey] = uxy(P(1,Node_P_index)+dx/2,P(2,Node_P_index),uParas);
            [u_sx,u_sy] = uxy(P(1,Node_P_index),P(2,Node_P_index)-dy/2,uParas);
            [u_nx,u_ny] = uxy(P(1,Node_P_index),P(2,Node_P_index)+dy/2,uParas);
            
            
            
            F_w=rho*u_wx*dy;
            F_e=rho*u_ex*dy;
            F_s=rho*u_sy*dx;
            F_n=rho*u_ny*dx;

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
    



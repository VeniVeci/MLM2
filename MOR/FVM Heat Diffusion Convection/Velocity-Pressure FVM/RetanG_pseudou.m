function [u_hat,Grid] = RetanG_pseudou(Grid)
%% retangle grid slover for X-Momentum Equation
% retangle grid slover for properties such as u, v, phi.
%
% Modifications:
% 14-May-2014, WeiX, first edition 

%%  
rho=Grid.rho;
k=Grid.k;
dx=Grid.dx;
dy=Grid.dy;
Num_x_u=Grid.Num_x_u;
Num_y_u=Grid.Num_y_u;
u_hat=Grid.u;

% u_P_index=0;
for x= 2: Num_x_u-1
    for y= 1: Num_y_u

%         u_P_index=u_P_index+1;            % point index 
%         u_W_index=u_P_index-Num_y_u;                  % identify west  neighbour (in a retangular grid)
%         u_E_index=u_P_index+Num_y_u;                  % identify east  neighbour (in a retangular grid)
%         u_S_index=u_P_index-1;                        % identify south neighbour (in a retangular grid)
%         u_N_index=u_P_index+1;                        % identify north neighbour (in a retangular grid)




        u_star_w=(Grid.u(x,y)+Grid.u(x-1,y))/2;
        u_star_e=(Grid.u(x,y)+Grid.u(x+1,y))/2;
        v_star_s=(Grid.v(x-1,y)+Grid.v(x,y))/2;
        v_star_n=(Grid.v(x-1,y+1)+Grid.v(x,y+1))/2;

        F_w= u_star_w*rho*dy;
        F_e= u_star_e*rho*dy;
        F_s= v_star_s*rho*dx;  
        F_n= v_star_n*rho*dx;  

        D_w= k*dy/dx;
        D_e= k*dy/dx;
        D_s= k*dx/dy;  
        D_n= k*dx/dy;       

        a_W=D_w+max([0, F_w]);                          % set default coefficients which is used for a inner node.
        a_E=D_e+max([0,-F_e]);
        a_S=D_s+max([0, F_s]);
        a_N=D_n+max([0,-F_n]);

        S_PW=0;
        S_PE=0;
        S_PS=0;
        S_PN=0;

        % Identify the boundary points------------------------------------
        if x==2                     % point on the left (next to boundary)
            u_star_W=0;             % ONLY for non slip condition
%             u_W_index=0;
        else
            u_star_W=Grid.u(x-1,y);
        end          

        if x==Num_x_u-1     % point on the right
            u_star_E=0;             % ONLY for non slip condition  
        else
            u_star_E=Grid.u(x+1,y);
        end         

        if y==1                     % point on the bottom   
            u_star_S=2*Grid.u_B-Grid.u(x,y);
        else
            u_star_S=Grid.u(x,y-1);            
        end         

        if y==Num_y_u     % point on the top 
            u_star_N=2*Grid.u_T-Grid.u(x,y);
        else
            u_star_N=Grid.u(x,y+1);         
            
        end        
        
        S_P=S_PW+S_PE+S_PS+S_PN;    
        a_P=S_P-a_W-a_E-a_S-a_N+F_w-F_e+F_s-F_n;
        b=0;
        
        u_hat(x,y)=(a_W*u_star_W+a_E*u_star_E+a_S*u_star_S+a_N*u_star_N+b)/(-a_P);
        Grid.a_P_u(x,y)=-a_P;           % save all the a_P information for later u-correction.
        
    end
end

end
function [v_hat,Grid] = RetanG_pseudov(Grid)
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
Num_x_v=Grid.Num_x_v;
Num_y_v=Grid.Num_y_v;
v_hat=Grid.v;


% v_P_index=0;
for x= 1: Num_x_v
    for y= 2: Num_y_v-1
        
%         v_P_index=v_P_index+1;                       % point index
%         v_W_index=v_P_index-(Num_y_v-2);   % identify west  neighbour (in a retangular grid)
%         v_E_index=v_P_index+(Num_y_v-2);   % identify east  neighbour (in a retangular grid)
%         v_S_index=v_P_index-1;                       % identify south neighbour (in a retangular grid)
%         v_N_index=v_P_index+1;                       % identify north neighbour (in a retangular grid)
        
%         p_star_s=Grid.p(x,y-1);
%         p_star_n=Grid.p(x,y);
        
        u_star_w=(Grid.u(x,y)+Grid.u(x,y-1))/2;
        u_star_e=(Grid.u(x+1,y)+Grid.u(x+1,y-1))/2;
        v_star_s=(Grid.v(x,y)+Grid.v(x,y-1))/2;
        v_star_n=(Grid.v(x,y)+Grid.v(x,y+1))/2;
              
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
        
%         a_W=D_w+F_w/2;                          % set default coefficients which is used for a inner node.
%         a_E=D_e-F_w/2;
%         a_S=D_s+F_s/2;
%         a_N=D_n-F_n/2;
        
        S_PW=0;
        S_PE=0;
        S_PS=0;
        S_PN=0;
        
        % Identify the boundary points------------------------------------
        if x==1                     % point on the left (next to boundary)
%            a_W=0;                  % set coefficient a_W zero
%            S_PW=-2*D_w-F_w;        % sources form the west boundary 
%            v_W_index=0;         % west neighbour is void 
             v_star_W=2*Grid.v_L-Grid.v(x,y);
        else
             v_star_W=Grid.v(x-1,y);
        end          
        
        if x==Num_x_v     % point on the right
%            a_E=0;
%            S_PE=-2*D_e+F_e; 
%            v_E_index=0;      
             v_star_E=2*Grid.v_R-Grid.v(x,y);
         else
            v_star_E=Grid.v(x+1,y);
        end         
        
        if y==2                     % point on the top
%             a_N=0;
%             S_PN=-2*D_n-F_n;
            v_star_S=0;  
        else
            v_star_S=Grid.v(x,y-1);
        end         
        
        if y==Num_y_v-1     % point on the bottom
%             a_S=0;
%             S_PS=-2*D_s+F_s;
            v_star_N=0;      
        else
            v_star_N=Grid.v(x,y+1);
        end       
        
        S_P=S_PW+S_PE+S_PS+S_PN;    
        a_P=S_P-a_W-a_E-a_S-a_N+F_w-F_e+F_s-F_n;
        b=0;
        
        v_hat(x,y)=(a_W*v_star_W+a_E*v_star_E+a_S*v_star_S+a_N*v_star_N+b)/(-a_P);
        Grid.a_P_v(x,y)=-a_P;           % save all the a_P information for later u-correction.
    end
end




end
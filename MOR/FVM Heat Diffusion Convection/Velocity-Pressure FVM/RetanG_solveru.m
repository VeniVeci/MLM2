function [u,Grid] = RetanG_solveru(Grid,option)
%% retangle grid slover for X-Momentum Equation
% retangle grid slover for properties such as u, v, phi.
%
% Modifications:
% 12-May-2014, WeiX, first edition 

%%  
%%  
% set default option
switch nargin
    case 1
       option='bkdiv';          %opetion is not specified. Using default setting  
end

rho=Grid.rho;
k=Grid.k;
dx=Grid.dx;
dy=Grid.dy;
Num_x_u=Grid.Num_x_u;
Num_y_u=Grid.Num_y_u;
u=Grid.u;

A=zeros(Num_y_u*(Num_x_u-2));     %left and right u are zero by default in incompressible fluid in no-slip wall container.
B=zeros(Num_y_u*(Num_x_u-2),1);
a_P_u=zeros(Num_x_u,Num_y_u);


u_P_index=0;
for x= 2: Num_x_u-1
    for y= 1: Num_y_u

        u_P_index=u_P_index+1;            % point index 
        u_W_index=u_P_index-Num_y_u;                  % identify west  neighbour (in a retangular grid)
        u_E_index=u_P_index+Num_y_u;                  % identify east  neighbour (in a retangular grid)
        u_S_index=u_P_index-1;                        % identify south neighbour (in a retangular grid)
        u_N_index=u_P_index+1;                        % identify north neighbour (in a retangular grid)

        p_star_w=Grid.p(x-1,y);
        p_star_e=Grid.p(x,y);

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
%           a_W=0;                  % set coefficient a_W zero
%           S_PW=-2*D_w-F_w;        % sources form the west boundary 
            u_W_index=0;            % west neighbour is void 
        end          

        if x==Num_x_u-1     % point on the right
%           a_E=0;
%           S_PE=-2*D_e+F_e; 
            u_E_index=0;      
        end         

        if y==1                     % point on the bottom
            a_S=0;
            S_PS=-2*D_s-F_s;
            u_S_index=0;      
        end         

        if y==Num_y_u     % point on the top
            a_N=0;
            S_PN=-2*D_n+F_n;
            u_N_index=0;      
        end        

        P_U=(p_star_w-p_star_e)*dy; 
        S_P=S_PW+S_PE+S_PS+S_PN;    
        S_U=-S_PW*Grid.u_L-S_PE*Grid.u_R-S_PS*Grid.u_B-S_PN*Grid.u_T;
        a_P=S_P-a_W-a_E-a_S-a_N+F_w-F_e+F_s-F_n;

        A(u_P_index,u_P_index)=a_P;        % putting coefficients in matrix A
        if u_W_index~=0   A(u_P_index,u_W_index)=a_W;end
        if u_E_index~=0   A(u_P_index,u_E_index)=a_E;end
        if u_S_index~=0   A(u_P_index,u_S_index)=a_S;end
        if u_N_index~=0   A(u_P_index,u_N_index)=a_N;end

        B(u_P_index,1)=-S_U-P_U;        % putting coefficients in matrix b
        Grid.a_P_u(x,y)=-a_P;           % save all the a_P information for later u-correction.
    end
end

switch option
    case 'bkdiv'
        u_temp=A\B;
    case 'TDMA'
        u_temp=TDMAsolver(A,B);
    case 'lsor'
        u_temp=lsqr(A,B);
    case 'penta'
        u_temp=pentsolve(A,B);
    case 'pcg'
        u_temp=pcg(A,B);
end

% u_temp=A\B;
% u_temp=TDMAsolver(A,B);


u_P_index=0;
for x= 2: Num_x_u-1
    for y= 1: Num_y_u
        u_P_index=u_P_index+1;            % point index 
        u(x,y)=u_temp(u_P_index,1);       
    end
end




end
function [p_correct] = RetanG_solverp(Grid,option)
%% retangle grid slover for corrected pressure
%
% Modifications:
% 12-May-2014, WeiX, first edition 

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
Num_x_p=Grid.Num_x_p;
Num_y_p=Grid.Num_y_p;
v=Grid.v;

A=zeros(Num_y_p*Num_x_p);
B=zeros(Num_y_p*Num_x_p,1);

p_P_index=0;
for x= 1: Num_x_p
    for y= 1: Num_y_p    
     
        p_P_index=p_P_index+1;                   % point index
        p_W_index=p_P_index-Num_y_p;             % identify west  neighbour (in a retangular grid)
        p_E_index=p_P_index+Num_y_p;             % identify east  neighbour (in a retangular grid)
        p_S_index=p_P_index-1;                   % identify south neighbour (in a retangular grid)
        p_N_index=p_P_index+1;                   % identify north neighbour (in a retangular gridp
        
        
        d_w=dy/Grid.a_P_u(x,y);
        d_e=dy/Grid.a_P_u(x+1,y);
        d_s=dx/Grid.a_P_v(x,y);
        d_n=dx/Grid.a_P_v(x,y+1);
        
        a_W=rho*dy*d_w;                       
        a_E=rho*dy*d_e;
        a_S=rho*dx*d_s;
        a_N=rho*dx*d_n;
        
        u_star_w=Grid.u(x,y);
        u_star_e=Grid.u(x+1,y);
        v_star_s=Grid.v(x,y);
        v_star_n=Grid.v(x,y+1);
        
        b_W=rho*dy*u_star_w;
        b_E=rho*dy*u_star_e;
        b_S=rho*dx*v_star_s;
        b_N=rho*dx*v_star_n;
         
        if x==1
            a_W=0;
            b_W=0;
            p_W_index=0;
        end
        
        if x==Num_x_p
            a_E=0;
            b_E=0;
            p_E_index=0;
        end
        
        if y==1          
            a_S=0;
            b_S=0;
            p_S_index=0;
        end
        
        if y==Num_y_p
            a_N=0;
            b_N=0;
            p_N_index=0;
        end
        
        a_P=-a_W-a_E-a_S-a_N;
        b=b_W-b_E+b_S-b_N;
        
        
        A(p_P_index,p_P_index)=a_P;        % putting coefficients in matrix A
        if p_W_index~=0   A(p_P_index,p_W_index)=a_W;end
        if p_E_index~=0   A(p_P_index,p_E_index)=a_E;end
        if p_S_index~=0   A(p_P_index,p_S_index)=a_S;end
        if p_N_index~=0   A(p_P_index,p_N_index)=a_N;end
        
        B(p_P_index,1)=-b;                 % putting coefficients in matrix b
        
    end
end

switch option
    case 'bkdiv'
        p_temp=A\B;
    case 'TDMA'
        p_temp=TDMAsolver(A,B);
    case 'lsor'
        p_temp=lsqr(A,B);
    case 'penta'
        p_temp=pentsolve(A,B);
    case 'pcg'
        p_temp=pcg(A,B);
end


p_correct=zeros(Num_x_p,Num_y_p);
p_P_index=0;
for x= 1: Num_x_p
    for y= 1: Num_y_p    
        p_P_index=p_P_index+1;      
        p_correct(x,y)=p_temp(p_P_index,1);

    end
end




end
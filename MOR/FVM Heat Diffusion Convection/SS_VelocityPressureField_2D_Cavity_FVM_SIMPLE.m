function SS_VelocityPressureField_2D_Cavity_FVM_SIMPLE
%Simulating the Velosity and Pressure field 2D equation  
...by the Finite Volume Method using SIMPLE algorithm (up-wind scheme)
...in a square plane in uniform grid
%    
% Governning equation : 
% 
% 
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
% 25-April-2014, WeiX, first edition 

clear
%% -----------------------Setup Parameters---------------------------------
%  -----------------------Problem Parameters---------------------------
Lx=1;                   % length of x axis is 1m
Ly=1;                   % length of y axis is 1m
rho=1000;                 % flow density
k=1;                    % ???????????

u_T=1;  v_T=0;            % Vertical & horizontal velocity of Top    lid                 
u_B=0;  v_B=0;            % Vertical & horizontal velocity of Bottom lid
u_L=0;  v_L=0;            % Vertical & horizontal velocity of Left   lid   
u_R=0;  v_R=0;            % Vertical & horizontal velocity of Right  lid   

% -----------------------Solver Parameters---------------------------
Scale_grid_x=50;        % number of volumes on x direction
Scale_grid_y=50;        % number of volumes on y direction

%% --------------Initializing computational Parameters------------------
dx=Lx/(Scale_grid_x-1);
dy=Ly/(Scale_grid_y-1);

Num_x_p=Scale_grid_x-1;     %number of the matrix of pressure of x
Num_y_p=Scale_grid_y-1;     %number of the matrix of pressure of y

Num_x_u=Scale_grid_x;
Num_y_u=Scale_grid_y-1;

Num_x_v=Scale_grid_x-1;
Num_y_v=Scale_grid_y;

% Grid=zeros(Scale_grid_y*2-1,Scale_grid_x*2-1);
p=zeros(Num_x_p,Num_y_p);
u=zeros(Num_x_u,Num_y_u);
u(1,:)=u_L;
u(end,:)=u_R;
v=zeros(Num_x_v,Num_y_v);
v(:,1)=v_B;
v(:,end)=v_T;

%------------------------------
p_star=rand(Num_x_p,Num_y_p);
u_star=rand(Num_x_u,Num_y_u);
u_star(1,:)=u_L;
u_star(end,:)=u_R;
v_star=rand(Num_x_v,Num_y_v);
v_star(:,1)=v_B;
v_star(:,end)=v_T;

%------------------------------
p_star=zeros(Num_x_p,Num_y_p)*10000;
u_star=zeros(Num_x_u,Num_y_u);
u_star(1,:)=u_L;
u_star(end,:)=u_R;
v_star=zeros(Num_x_v,Num_y_v);
v_star(:,1)=v_B;
v_star(:,end)=v_T;


%% -----------------------F V M Solver-------------------------------------
p_record=[];
u_record=[];
for i=1:50
%% ---------------------Solving X-Momentum Equation------------------------
A=zeros(Num_y_u*(Num_x_u-2));
B=zeros(Num_y_u*(Num_x_u-2),1);
a_P_u=zeros(Num_x_u,Num_y_u);
u_P_index=0;
for x= 2: Num_x_u-1
    for y= 1: Num_y_u
        
        u_P_index=u_P_index+1;            % point index 
        
        % Information Record not useful yet 
%         u_location(1,Node_u_index)=dx*(x-1);             % u_location(1,:) is the X cooridnate of points
%         u_location(2,Node_u_index)=dy*(1/2+(y-1));       % u_location(2,:) is the y cooridnate of points
%         u_location(3,Node_u_index)=Node_u_index;         % u_location(3,:) is the index of points
                     
        u_W_index=u_P_index-Num_y_u;                  % identify west  neighbour (in a retangular grid)
        u_E_index=u_P_index+Num_y_u;                  % identify east  neighbour (in a retangular grid)
        u_S_index=u_P_index-1;                        % identify south neighbour (in a retangular grid)
        u_N_index=u_P_index+1;                        % identify north neighbour (in a retangular grid)
        
    
        p_star_w=p_star(x-1,y);
        p_star_e=p_star(x,y);
       
        u_star_w=(u_star(x,y)+u_star(x-1,y))/2;
        u_star_e=(u_star(x,y)+u_star(x+1,y))/2;
        v_star_s=(v_star(x-1,y)+v_star(x,y))/2;
        v_star_n=(v_star(x-1,y+1)+v_star(x,y+1))/2;
        
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
        S_U=-S_PW*u_L-S_PE*u_R-S_PS*u_B-S_PN*u_T;
        a_P=S_P-a_W-a_E-a_S-a_N+F_w-F_e+F_s-F_n;

        A(u_P_index,u_P_index)=a_P;        % putting coefficients in matrix A
        if u_W_index~=0   A(u_P_index,u_W_index)=a_W;end
        if u_E_index~=0   A(u_P_index,u_E_index)=a_E;end
        if u_S_index~=0   A(u_P_index,u_S_index)=a_S;end
        if u_N_index~=0   A(u_P_index,u_N_index)=a_N;end
        
        B(u_P_index,1)=-S_U-P_U;        % putting coefficients in matrix b
        a_P_u(x,y)=-a_P;                 % save all the a_P information for later u-correction.
        

    end
end

u_temp=A\B;
% u_temp=TDMAsolver(A,B);

u_P_index=0;
for x= 2: Num_x_u-1
    for y= 1: Num_y_u
        u_P_index=u_P_index+1;            % point index 
        u(x,y)=u_temp(u_P_index,1);       
    end
end


%% ---------------------Solving Y-Momentum Equation------------------------        
A=zeros(Num_x_v*(Num_y_v-2));
B=zeros(Num_x_v*(Num_y_v-2),1);
a_P_v=zeros(Num_x_v,Num_y_v);
v_P_index=0;

for x= 1: Num_x_v
    for y= 2: Num_y_v-1
        
        v_P_index=v_P_index+1;                       % point index
        v_W_index=v_P_index-(Num_y_v-2);   % identify west  neighbour (in a retangular grid)
        v_E_index=v_P_index+(Num_y_v-2);   % identify east  neighbour (in a retangular grid)
        v_S_index=v_P_index-1;                       % identify south neighbour (in a retangular grid)
        v_N_index=v_P_index+1;                       % identify north neighbour (in a retangular grid)
        
        p_star_s=p_star(x,y-1);
        p_star_n=p_star(x,y);
        
        u_star_w=(u_star(x,y)+u_star(x,y-1))/2;
        u_star_e=(u_star(x+1,y)+u_star(x+1,y-1))/2;
        v_star_s=(v_star(x,y)+v_star(x,y-1))/2;
        v_star_n=(v_star(x,y)+v_star(x,y+1))/2;
              
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
           a_W=0;                  % set coefficient a_W zero
           S_PW=-2*D_w-F_w;        % sources form the west boundary 
           v_W_index=0;         % west neighbour is void 
        end          
        
        if x==Num_x_v     % point on the right
           a_E=0;
           S_PE=-2*D_e+F_e; 
           v_E_index=0;      
        end         
        
        if y==2                     % point on the top
%             a_N=0;
%             S_PN=-2*D_n-F_n;
            v_S_index=0;      
        end         
        
        if y==Num_y_v-1     % point on the bottom
%             a_S=0;
%             S_PS=-2*D_s+F_s;
            v_N_index=0;      
        end       
        
        % Assembly of the linear equation---------------------------------
        P_U=(p_star_s-p_star_n)*dx; 
        S_P=S_PW+S_PE+S_PS+S_PN;    
        S_U=-S_PW*v_L-S_PE*v_R-S_PS*v_B-S_PN*v_T;
        a_P=S_P-a_W-a_E-a_S-a_N+F_w-F_e+F_s-F_n;

        A(v_P_index,v_P_index)=a_P;        % putting coefficients in matrix A
        if v_W_index~=0   A(v_P_index,v_W_index)=a_W;end
        if v_E_index~=0   A(v_P_index,v_E_index)=a_E;end
        if v_S_index~=0   A(v_P_index,v_S_index)=a_S;end
        if v_N_index~=0   A(v_P_index,v_N_index)=a_N;end
           
        B(v_P_index,1)=-S_U-P_U;        % putting coefficients in matrix b
        
        a_P_v(x,y)=-a_P;

    end
end

v_temp=A\B;
% v_temp=TDMAsolver(A,B);

v_P_index=0;
for x= 1: Num_x_v
    for y= 2: Num_y_v-1
        v_P_index=v_P_index+1;            % point index 
        v(x,y)=v_temp(v_P_index,1);       
    end
end

%% ----------------Update uv field-----------------------------------------
u_star=u;
v_star=v;

%% ----------------Solving Pressure-Correction Equation--------------------

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
        
        
        d_w=dy/a_P_u(x,y);
        d_e=dy/a_P_u(x+1,y);
        d_s=dx/a_P_v(x,y);
        d_n=dx/a_P_v(x,y+1);
        
        a_W=rho*dy*d_w;                       
        a_E=rho*dy*d_e;
        a_S=rho*dx*d_s;
        a_N=rho*dx*d_n;
        
        u_star_w=u_star(x,y);
        u_star_e=u_star(x+1,y);
        v_star_s=v_star(x,y);
        v_star_n=v_star(x,y+1);
        
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

% p_temp = bicg(A,B);                 %fail
% p_temp = pcg(A,B);                 %fail
% p_temp = lsqr(A,B);                 % work
% p_temp = cgs(A,B);                %fail
% p_temp = gmres(A,B);              %fail
p_temp=TDMAsolver(A,B);           % work
% p_temp=LESsolver(A,B,'BSOR');     %fail
norm(A*p_temp-B)
% p_temp = linsolve(A,B);           %fail

% p_temp=A\B;
% p_temp=inv(A)*B;
% det(A)
% p_correct=reshape(p_temp,[Grid.Info.Num_y_p,Grid.Info.Num_x_p]);

p_correct=zeros(Num_x_p,Num_y_p);
p_P_index=0;
for x= 1: Num_x_p
    for y= 1: Num_y_p    
        p_P_index=p_P_index+1;
        
%         if p_temp(p_P_index,1) < 0 p_temp(p_P_index,1)=0;end
            
        p_correct(x,y)=p_temp(p_P_index,1);

    end
end

%% ------------------Correct uvp fields------------------------------------
Alpha_u=1;
Alph1_v=1;
Alph1_p=0.1;
% ----------p field-------------------------------
p_star=p_star+p_correct*Alph1_p;

%check if the p filed is correct

% %---method 1------------
% minip=min(p_star(:));
% p_star=p_star-minip*ones(Num_x_p,Num_y_p);

%---method 2------------
% for x= 1: Num_x_p
%     for y= 1: Num_y_p    
%      if p_star(x,y) < 0 p_star(x,y)=0;end
%     end
% end

% ----------u field-------------------------------
for x= 2: Num_x_u-1
    for y= 1: Num_y_u
        p_star_w=p_correct(x-1,y);
        p_star_e=p_correct(x,y);
        u_star(x,y)=u_star(x,y)+(p_star_w-p_star_e)*dy/a_P_u(x,y);     
    end
end
% ----------v field-------------------------------
for x= 1: Num_x_v
    for y= 2: Num_y_v-1
        
        p_star_s=p_correct(x,y-1);
        p_star_n=p_correct(x,y);
        v_star(x,y)=v_star(x,y)+(p_star_s-p_star_n)*dx/a_P_v(x,y);     

    end
end

%% --------Convergence Check----------------------------------------------
p_abs=sum(abs(p_star(:)));
p_record=[p_record;p_abs];

u_abs=sum(abs(u_star(:)));
u_record=[u_record;u_abs];


end

%% ------------Visualization----------------------------------------------

 uplot(1:Num_x_p,1:Num_y_p)=0.5*(u(1:Num_x_u-1,1:Num_y_u)+u(2:Num_x_u,1:Num_y_u));
 vplot(1:Num_x_p,1:Num_y_p)=0.5*(v(1:Num_x_v,1:Num_y_v-1)+v(1:Num_x_v,2:Num_y_v));
for x= 1:Num_x_p
    for y=1:Num_y_p
        corx(x,y)=x;
        cory(x,y)=y;
    end
end

 quiver(corx,cory,uplot,vplot,1,'k-');

end


function X=TDMAsolver(A,b)

m=length(b);                 % m is the number of rows
X=zeros(m,1);
A(1,2)= A(1,2)  ./ A(1,1);    % Division by zero risk.
b(1)=  b(1)    ./ A(1,1);    % Division by zero would imply a singular matrix

for i=2:m-1
    temp=  A(i,i) - A(i,i-1) .* A(i-1,i);
    A(i,i+1)=  A(i,i+1)  ./ temp;
    b(i)= ( b(i) - A(i,i-1) .* b(i-1) )  ./ temp;
end 

i=m;
X(m)=(b(i) - A(i,i-1) .* b(i-1))  ./ (A(i,i) - A(i,i-1) .* A(i-1,i));

for i=m-1:-1:1
X(i)=  -A(i,i+1) .* X(i+1) + b(i);
end
end



function SS_VelocityPressureField_2D_Cavity_FVM_SIMPLE_MG
%Simulating the Velosity and Pressure field 2D equation  
...by the Finite Volume Method using SIMPLE algorithm (up-wind scheme)
...in a square plane in uniform grid in a "merge grid" n
%    
% Governning equation : 
% 
% 
% The point are discritized as follow grid. e.g a 3*4 grid. Numbers are
% index of point
% 
%        1--5--9
%        |  |  | 
%        2--6--10
%        |  |  |
%     ^  3--7--11
%     |  |  |  |
%  u_y|  4--8--12
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
rho=1;                 % flow density
k=0.1;                    % ?????????????????????????????????????

u_T=1;  v_T=0;            % Vertical & horizontal velocity of Top    lid                 
u_B=0;  v_B=0;            % Vertical & horizontal velocity of Bottom lid
u_L=0;  v_L=0;            % Vertical & horizontal velocity of Left   lid   
u_R=0;  v_R=0;            % Vertical & horizontal velocity of Right  lid  

% -----------------------Solver Parameters---------------------------
Scale_grid_x=20;        % number of volumes on x direction
Scale_grid_y=20;        % number of volumes on y direction

%% --------------Initializing Mesh----------------------------------------
[Grid] = UniG_Maker(Scale_grid_x,Scale_grid_y,Lx,Ly);

%% --------------Guess uvp fields------------------------------------------
u_star=rand(Grid.Info.Num_y_u,Grid.Info.Num_x_u);
u_star(:,1)=u_L;
u_star(:,end)=u_R;

v_star=rand(Grid.Info.Num_y_v,Grid.Info.Num_x_v);
v_star(end,:)=v_B;
v_star(1,:)=v_T;

p_star=rand(Grid.Info.Num_y_p,Grid.Info.Num_x_p);

% [Grid] = UniG_Update(Grid,u_star,v_star,p_star);
[Grid] = UniG_UpdateSolo(Grid,u_star,'u');
[Grid] = UniG_UpdateSolo(Grid,v_star,'v');
[Grid] = UniG_UpdateSolo(Grid,p_star,'p');
%% --------------Test area-------------------------------------------------
u_star=ones(Grid.Info.Num_y_u,Grid.Info.Num_x_u);
u_star(:,1)=u_L;
u_star(:,end)=u_R;

v_star=ones(Grid.Info.Num_y_v,Grid.Info.Num_x_v);
v_star(end,:)=v_B;
v_star(1,:)=v_T;

p_star=ones(Grid.Info.Num_y_p,Grid.Info.Num_x_p);

[Grid] = UniG_UpdateSolo(Grid,u_star,'u');
[Grid] = UniG_UpdateSolo(Grid,v_star,'v');
[Grid] = UniG_UpdateSolo(Grid,p_star,'p');

%% -----------------------F V M Solver-------------------------------------
%%--test code
% [Grid1] = UniG_Maker(5,6,Lx,Ly);
% 
% u1=1:25;
% u1 = reshape(u1,[5,5])
% 
% v1=100:100:2400;
% v1 = reshape(v1,[6,4])
% 
% [Grid1] = UniG_Update(Grid1,u1,v1)

p_record=[];
u_record=[];
v_record=[];

for i=1:100

%% ---------------------Solving X-Momentum Equation------------------------
dx=Grid.Info.dx;
dy=Grid.Info.dy;
A=zeros(Grid.Info.Num_y_u*(Grid.Info.Num_x_u-2));
B=zeros(Grid.Info.Num_y_u*(Grid.Info.Num_x_u-2),1);
a_P_u=zeros(Grid.Info.Num_y_u,Grid.Info.Num_x_u);
a_P_Gu=zeros(Grid.Info.Num_y_full,Grid.Info.Num_x_full);

u_P_index=0;
for x= 2: Grid.Info.Num_x_u-1
    for y= 1: Grid.Info.Num_y_u    
        
        u_P_index=u_P_index+1;                   % point index
        u_W_index=u_P_index-Grid.Info.Num_y_u;   % identify west  neighbour (in a retangular grid)
        u_E_index=u_P_index+Grid.Info.Num_y_u;   % identify east  neighbour (in a retangular grid)
        u_S_index=u_P_index+1;                   % identify south neighbour (in a retangular grid)
        u_N_index=u_P_index-1;                   % identify north neighbour (in a retangular grid)
        
        [G_x,G_y] = uvp2GP_UniG(x,y,'u');        % identify global coordnate G_x,G_y
        
        p_star_w=Grid.globalp(G_y,G_x-1);
        p_star_e=Grid.globalp(G_y,G_x+1);
        
        u_star_w=Grid.globalu(G_y,G_x-1);
        u_star_e=Grid.globalu(G_y,G_x+1);
        v_star_s=Grid.globalv(G_y+1,G_x);
        v_star_n=Grid.globalv(G_y-1,G_x);

        F_w= u_star_w*rho*dy;
        F_e= u_star_e*rho*dy;
        F_s= v_star_s*rho*dx;  
        F_n= v_star_n*rho*dx;  
        
        D_w= k*dy/dx;
        D_e= k*dy/dx;
        D_s= k*dx/dy;  
        D_n= k*dx/dy;  
       
        Area_p=dy;
        
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
            u_W_index=0;         % west neighbour is void 
        end          
        
        if x==Grid.Info.Num_x_u-1     % point on the right
%           a_E=0;
%           S_PE=-2*D_e+F_e; 
            u_E_index=0;      
        end         
        
        if y==1                     % point on the top
            a_N=0;
            S_PN=-2*D_n+F_n;
            u_N_index=0;      
        end         
        
        if y==Grid.Info.Num_y_u     % point on the bottom
            a_S=0;
            S_PS=-2*D_s-F_s;
            u_S_index=0;      
        end        
          
        % Assembly of the linear equation---------------------------------
        P_U=(p_star_w-p_star_e)*Area_p; 
        S_P=S_PW+S_PE+S_PS+S_PN;    
        S_U=-S_PW*u_L-S_PE*u_R-S_PS*u_B-S_PN*u_T;
        a_P=S_P-a_W-a_E-a_S-a_N+F_w-F_e+F_s-F_n;

        A(u_P_index,u_P_index)=a_P;        % putting coefficients in matrix A
        if u_W_index~=0   A(u_P_index,u_W_index)=a_W;end
        if u_E_index~=0   A(u_P_index,u_E_index)=a_E;end
        if u_S_index~=0   A(u_P_index,u_S_index)=a_S;end
        if u_N_index~=0   A(u_P_index,u_N_index)=a_N;end
           
        B(u_P_index,1)=-S_U-P_U;        % putting coefficients in matrix b
        
        a_P_u(y,x)=a_P;                 % save all the a_P information for later u-correction.
        a_P_Gu(G_y,G_x)=-a_P;
    end
end
u_temp=A\B;
% u_temp=inv(A)*B;
% u_temp=TDMAsolver(A,B);

u_Inner=reshape(u_temp,[Grid.Info.Num_y_u,Grid.Info.Num_x_u-2]);
u_new=[zeros(Grid.Info.Num_y_u,1),u_Inner,zeros(Grid.Info.Num_y_u,1)];

%% ---------------------Solving Y-Momentum Equation------------------------
dx=Grid.Info.dx;
dy=Grid.Info.dy;
A=zeros(Grid.Info.Num_x_v*(Grid.Info.Num_y_v-2));
B=zeros(Grid.Info.Num_x_v*(Grid.Info.Num_y_v-2),1);
a_P_v=zeros(Grid.Info.Num_y_v,Grid.Info.Num_x_v-2);
a_P_Gv=zeros(Grid.Info.Num_y_full,Grid.Info.Num_x_full);

v_P_index=0;
for x= 1: Grid.Info.Num_x_v
    for y= 2: Grid.Info.Num_y_v-1    
        
        v_P_index=v_P_index+1;                   % point index
        v_W_index=v_P_index-(Grid.Info.Num_y_v-2);   % identify west  neighbour (in a retangular grid)
        v_E_index=v_P_index+(Grid.Info.Num_y_v-2);   % identify east  neighbour (in a retangular grid)
        v_S_index=v_P_index+1;                   % identify south neighbour (in a retangular grid)
        v_N_index=v_P_index-1;                   % identify north neighbour (in a retangular grid)
        
        [G_x,G_y] = uvp2GP_UniG(x,y,'v');        % identify global coordnate G_x,G_y
        
        p_star_s=Grid.globalp(G_y+1,G_x);
        p_star_n=Grid.globalp(G_y-1,G_x);
        
        u_star_w=Grid.globalu(G_y,G_x-1);
        u_star_e=Grid.globalu(G_y,G_x+1);
        v_star_s=Grid.globalv(G_y+1,G_x);
        v_star_n=Grid.globalv(G_y-1,G_x);

        F_w= u_star_w*rho*dy;
        F_e= u_star_e*rho*dy;
        F_s= v_star_s*rho*dx;  
        F_n= v_star_n*rho*dx;  
        
        D_w= k*dy/dx;
        D_e= k*dy/dx;
        D_s= k*dx/dy;  
        D_n= k*dx/dy;  
       
        Area_p=dx;
        
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
        
        if x==Grid.Info.Num_x_v     % point on the right
           a_E=0;
           S_PE=-2*D_e+F_e; 
           v_E_index=0;      
        end         
        
        if y==2                     % point on the top
%             a_N=0;
%             S_PN=-2*D_n-F_n;
            v_N_index=0;      
        end         
        
        if y==Grid.Info.Num_y_v-1     % point on the bottom
%             a_S=0;
%             S_PS=-2*D_s+F_s;
            v_S_index=0;      
        end        
          
        % Assembly of the linear equation---------------------------------
        P_U=(p_star_s-p_star_n)*Area_p; 
        S_P=S_PW+S_PE+S_PS+S_PN;    
        S_U=-S_PW*v_L-S_PE*v_R-S_PS*v_B-S_PN*v_T;
        a_P=S_P-a_W-a_E-a_S-a_N+F_w-F_e+F_s-F_n;

        A(v_P_index,v_P_index)=a_P;        % putting coefficients in matrix A
        if v_W_index~=0   A(v_P_index,v_W_index)=a_W;end
        if v_E_index~=0   A(v_P_index,v_E_index)=a_E;end
        if v_S_index~=0   A(v_P_index,v_S_index)=a_S;end
        if v_N_index~=0   A(v_P_index,v_N_index)=a_N;end
           
        B(v_P_index,1)=-S_U-P_U;        % putting coefficients in matrix b
        
        a_P_v(y,x)=a_P;
        a_P_Gv(G_y,G_x)=-a_P;
    end
end
v_temp=A\B;
% p_temp=inv(A)*B;
% v_temp=TDMAsolver(A,B);

v_Inner=reshape(v_temp,[Grid.Info.Num_y_v-2,Grid.Info.Num_x_v]);
v_new=[zeros(1,Grid.Info.Num_x_v);v_Inner;zeros(1,Grid.Info.Num_x_v)];

%% ----------------Update Velocity Field-----------------------------------
[Grid] = UniG_UpdateSolo(Grid,u_new,'u');
[Grid] = UniG_UpdateSolo(Grid,v_new,'v');

%% ----------------Solving Pressure-Correction Equation--------------------
dx=Grid.Info.dx;
dy=Grid.Info.dy;
A=zeros(Grid.Info.Num_y_p*Grid.Info.Num_x_p);
B=zeros(Grid.Info.Num_y_p*Grid.Info.Num_x_p,1);

p_P_index=0;
for x= 1: Grid.Info.Num_x_p
    for y= 1: Grid.Info.Num_y_p    
     
        p_P_index=p_P_index+1;                   % point index
        p_W_index=p_P_index-Grid.Info.Num_y_p;   % identify west  neighbour (in a retangular grid)
        p_E_index=p_P_index+Grid.Info.Num_y_p;   % identify east  neighbour (in a retangular grid)
        p_S_index=p_P_index+1;                   % identify south neighbour (in a retangular grid)
        p_N_index=p_P_index-1;                   % identify north neighbour (in a retangular gridp
        
        [G_x,G_y] = uvp2GP_UniG(x,y,'p');        % identify global coordnate G_x,G_y    
           
%         u_star_w=Grid.globalu(G_y,G_x-1);
%         u_star_e=Grid.globalu(G_y,G_x+1);
%         v_star_s=Grid.globalv(G_y+1,G_x);
%         v_star_n=Grid.globalv(G_y-1,G_x);
        
        d_w=dy/a_P_Gu(G_y,G_x-1);
        d_e=dy/a_P_Gu(G_y,G_x+1);
        d_s=dx/a_P_Gv(G_y+1,G_x);
        d_n=dx/a_P_Gv(G_y-1,G_x);
        
        a_W=rho*dy*d_w;                       
        a_E=rho*dy*d_e;
        a_S=rho*dx*d_s;
        a_N=rho*dx*d_n;
        
        u_star_w=Grid.globalu(G_y,G_x-1);
        u_star_e=Grid.globalu(G_y,G_x+1);
        v_star_s=Grid.globalv(G_y+1,G_x);
        v_star_n=Grid.globalv(G_y-1,G_x);
        
        b_W=rho*dy*u_star_w;
        b_E=rho*dy*u_star_e;
        b_S=rho*dx*v_star_s;
        b_N=rho*dx*v_star_n;
         
        if x==1
            a_W=0;
            b_W=0;
            p_W_index=0;
        end
        
        if x==Grid.Info.Num_x_p
            a_E=0;
            b_E=0;
            p_E_index=0;
        end
        
        if y==1
            a_N=0;
            b_N=0;
            p_N_index=0;
        end
        
        if y==Grid.Info.Num_y_p
            a_S=0;
            b_S=0;
            p_S_index=0;
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

% p_temp=inv(A)*B;
% p_temp=A\B;
p_temp=TDMAsolver(A,B);

% det(A)
p_correct=reshape(p_temp,[Grid.Info.Num_y_p,Grid.Info.Num_x_p]);

%% ------------------Correct uvp fields------------------------------------
% ----------p field-------------------------------
p_star=Grid.p+p_correct*0.1;
[Grid] = UniG_UpdateSolo(Grid,p_star,'p');

% expand p_correct to p_global_correct
p_G_correct=zeros(Grid.Info.Num_y_full,Grid.Info.Num_x_full);
for x=2:2:Grid.Info.Num_x_full-1
    for y = 2:2:Grid.Info.Num_y_full-1
        p_G_correct(y,x)=p_correct(y/2,x/2);
    end
end

% ----------u field-------------------------------

u_star=zeros(Grid.Info.Num_y_u,Grid.Info.Num_x_u);
for x= 2: Grid.Info.Num_x_u-1
    for y= 1: Grid.Info.Num_y_u    
        
        [G_x,G_y] = uvp2GP_UniG(x,y,'u');        % identify global coordnate G_x,G_y
        
        p_correct_w=p_G_correct(G_y,G_x-1);
        p_correct_e=p_G_correct(G_y,G_x+1);
%       Area=dy;       

        d=dy/a_P_Gu(G_y,G_x);
        u_star(y,x)=Grid.u(y,x)+d*(p_correct_w-p_correct_e);
        
    end
end

u_star(:,1)=0;
u_star(:,end)=0;
[Grid] = UniG_UpdateSolo(Grid,u_star,'u');

% ----------v field-------------------------------
v_star=zeros(Grid.Info.Num_y_v,Grid.Info.Num_x_v);
for x= 1: Grid.Info.Num_x_v
    for y= 2: Grid.Info.Num_y_v-1    
      
        [G_x,G_y] = uvp2GP_UniG(x,y,'v');        % identify global coordnate G_x,G_y
        
        p_correct_s=p_G_correct(G_y+1,G_x);
        p_correct_n=p_G_correct(G_y-1,G_x);
%       Area=dx;       
        d=dx/a_P_Gv(G_y,G_x);
%         d=dx/a_P_v(y,x);
        v_star(y,x)=Grid.v(y,x)+d*(p_correct_s-p_correct_n);
        
    end
end

v_star(1,:)=0;
v_star(end,:)=0;
[Grid] = UniG_UpdateSolo(Grid,v_star,'v');

%% --------Convergence Check----------------------------------------------
p_abs=sum(abs(Grid.p(:)));
p_record=[p_record;p_abs];

u_abs=sum(abs(u_star(:)));
u_record=[u_record;u_abs];

v_abs=sum(abs(v_star(:)));
v_record=[v_record;v_abs];

end

%% -----------------------------visualization------------------------------

% quiver(1:Grid.Info.Num_x_full,1:Grid.Info.Num_y_full,Grid.globalu,Grid.globalv,0.6,'k-');
uplot(1:Grid.Info.Num_y_p,1:Grid.Info.Num_x_p)=0.5*(u_star(1:Grid.Info.Num_y_u,1:Grid.Info.Num_x_u-1)+u_star(1:Grid.Info.Num_y_u,2:Grid.Info.Num_x_u));
vplot(1:Grid.Info.Num_y_p,1:Grid.Info.Num_x_p)=0.5*(v_star(1:Grid.Info.Num_y_v-1,1:Grid.Info.Num_x_v)+v_star(2:Grid.Info.Num_y_v,1:Grid.Info.Num_x_v));
for x= 1:Grid.Info.Num_x_p
    for y=1:Grid.Info.Num_y_p
        corx(y,x)=x;
        cory(y,x)=y;
    end
end
% quiver(corx,cory,uplot,vplot,0.6,'k-');
quiver(1:Grid.Info.Num_x_p,1:Grid.Info.Num_y_p,uplot,-vplot,1,'k-');


end

%% --------Functions-------------------------------------------------------

function [Grid] = UniG_Maker(grid_Num_x,grid_Num_y,Lx,Ly)
dx=Lx/(grid_Num_x-1);
dy=Ly/(grid_Num_y-1);

Num_x_p=grid_Num_x-1;     %number of the matrix of pressure of x
Num_y_p=grid_Num_y-1;     %number of the matrix of pressure of y

Num_x_u=grid_Num_x;
Num_y_u=grid_Num_y-1;

Num_x_v=grid_Num_x-1;
Num_y_v=grid_Num_y;

Grid.ponit=zeros(2*grid_Num_y-1,2*grid_Num_x-1);
Grid.globalu=zeros(2*grid_Num_y-1,2*grid_Num_x-1);
Grid.globalv=zeros(2*grid_Num_y-1,2*grid_Num_x-1);
Grid.globalp=zeros(2*grid_Num_y-1,2*grid_Num_x-1);
Grid.u=zeros(Num_y_u,Num_x_u);
Grid.v=zeros(Num_y_v,Num_x_v);
Grid.p=zeros(Num_y_p,Num_x_p);

Grid.Info.Num_x_full=2*grid_Num_x-1;
Grid.Info.Num_y_full=2*grid_Num_y-1;
Grid.Info.dx=dx;
Grid.Info.dy=dy;
Grid.Info.Num_x_p=Num_x_p;
Grid.Info.Num_y_p=Num_y_p;
Grid.Info.Num_x_u=Num_x_u;
Grid.Info.Num_y_u=Num_y_u;
Grid.Info.Num_x_v=Num_x_v;
Grid.Info.Num_y_v=Num_y_v;
end






function [Grid] = UniG_UpdateSolo(Grid,x,type)

switch type
    case 'u'
        u=x;
        [rowu,colu]=size(u);
        [rowGu,colGu]=size(Grid.u);
        
        if rowu==rowGu & colu==colGu
            Grid.u=u;
            %-------interpret glovalu field using linear combination of u--------------
            for x=1:2:Grid.Info.Num_x_full
                for y = 2:2:Grid.Info.Num_y_full-1
                    Grid.globalu(y,x)=u(y/2,(x+1)/2);
                end
            end
            for x=1:2:Grid.Info.Num_x_full
                for y = 3:2:Grid.Info.Num_y_full-2
                    Grid.globalu(y,x)=(Grid.globalu(y-1,x)+Grid.globalu(y+1,x))/2;
                end
            end
            for x=2:2:Grid.Info.Num_x_full-1
                for y = 2:2:Grid.Info.Num_y_full-1
                    Grid.globalu(y,x)=(Grid.globalu(y,x-1)+Grid.globalu(y,x+1))/2;
                end
            end
        else
            error('size of u doed not match')
        end
    case 'v'
        v=x;
        [rowv,colv]=size(v);
        [rowGv,colGv]=size(Grid.v);   
        
        if rowv==rowGv & colv==colGv
            
            Grid.v=v;
            %-------interpret glovalv field using linear combination of v--------------
            for x=2:2:Grid.Info.Num_x_full-1
                for y = 1:2:Grid.Info.Num_y_full
                    Grid.globalv(y,x)=v((y+1)/2,x/2);
                end
            end
            for x=2:2:Grid.Info.Num_x_full-1
                for y = 2:2:Grid.Info.Num_y_full-1
                    Grid.globalv(y,x)=(Grid.globalv(y-1,x)+Grid.globalv(y+1,x))/2;
                end
            end
            for x=3:2:Grid.Info.Num_x_full-2
                for y = 1:2:Grid.Info.Num_y_full
                    Grid.globalv(y,x)=(Grid.globalv(y,x-1)+Grid.globalv(y,x+1))/2;
                end
            end
            
        else
            error('size of v doed not match')
        end

    case 'p'
        p=x;
        [rowp,colp]=size(p);
        [rowGp,colGp]=size(Grid.p);   
        
        if rowp==rowGp & colp==colGp
            
            Grid.p=p;       
            %-------update pressure field (p) in a global grid -----------------------
            for x=2:2:Grid.Info.Num_x_full-1
                for y = 2:2:Grid.Info.Num_y_full-1
                    Grid.globalp(y,x)=p(y/2,x/2);
                end
            end
            
        else
            error('size of p doed not match')
        end
        
    otherwise
        error('No such type in UniGrid');
end

end


%------------------------------------------------------------------------


function [Grid] = UniG_Update(Grid,u,v,p)
Grid.u=u;
Grid.v=v;
Grid.p=p;
%-------interpret glovalu field using linear combination of u--------------
for x=1:2:Grid.Info.Num_x_full
    for y = 2:2:Grid.Info.Num_y_full-1
        Grid.globalu(y,x)=u(y/2,(x+1)/2);
    end
end
for x=1:2:Grid.Info.Num_x_full
    for y = 3:2:Grid.Info.Num_y_full-2
        Grid.globalu(y,x)=(Grid.globalu(y-1,x)+Grid.globalu(y+1,x))/2;
    end
end
for x=2:2:Grid.Info.Num_x_full-1
    for y = 2:2:Grid.Info.Num_y_full-1
        Grid.globalu(y,x)=(Grid.globalu(y,x-1)+Grid.globalu(y,x+1))/2;
    end
end

%-------interpret glovalv field using linear combination of v--------------
for x=2:2:Grid.Info.Num_x_full-1
    for y = 1:2:Grid.Info.Num_y_full
        Grid.globalv(y,x)=v((y+1)/2,x/2);
    end
end
for x=2:2:Grid.Info.Num_x_full-1
    for y = 2:2:Grid.Info.Num_y_full-1
        Grid.globalv(y,x)=(Grid.globalv(y-1,x)+Grid.globalv(y+1,x))/2;
    end
end
for x=3:2:Grid.Info.Num_x_full-2
    for y = 1:2:Grid.Info.Num_y_full
        Grid.globalv(y,x)=(Grid.globalv(y,x-1)+Grid.globalv(y,x+1))/2;
    end
end

%-------update pressure field (p) in a global grid -----------------------
for x=2:2:Grid.Info.Num_x_full-1
    for y = 2:2:Grid.Info.Num_y_full-1
        Grid.globalp(y,x)=p(y/2,x/2);
    end
end

end



function [x,y,Fullindex,SoluIndex,type] = GP2uvp_UniG(G_x,G_y,grid_Num_x,grid_Num_y)
%Transfer a global point location to a sepcific property(u, v, or p)
... in a uniform grid.
%SoluIndex are index of unknown, sulution requied point.
%GgINdex   are the global grid index

Num_x_p=grid_Num_x-1;     %number of the matrix of pressure of x
Num_y_p=grid_Num_y-1;     %number of the matrix of pressure of y

Num_x_u=grid_Num_x;
Num_y_u=grid_Num_y-1;

Num_x_v=grid_Num_x-1;
Num_y_v=grid_Num_y;


if (mod(G_x,2) == 0) &&  (mod(G_y,2) == 0)
    type='p';
    x=G_x/2;
    y=G_y/2;
    Fullindex=(x-1)*Num_y_p+y;
    SoluIndex=(x-1)*Num_y_p+y;
end

if (mod(G_x,2) == 0) &&  (mod(G_y,2) ~= 0)
    type='v';
    x=G_x/2;
    y=(G_y+1)/2;
    Fullindex=(x-1)*Num_y_v+y;
    SoluIndex=(x-1)*(Num_y_v-2)+y-1;
end

if (mod(G_x,2) ~= 0) &&  (mod(G_y,2) == 0)
    type='u';
    x=(G_x+1)/2;
    y=G_y/2;
    Fullindex=(x-1)*Num_y_u+y;
    SoluIndex=(x-2)*Num_y_u+y;
end

if (mod(G_x,2) ~= 0) &&  (mod(G_y,2) ~= 0)
    type='None';
    x=(G_x+1)/2;
    y=(G_y+1)/2;
    Fullindex=(x-1)*Num_y_u+y;
    SoluIndex=(x-1)*Num_y_u+y;
end
end



function [G_x,G_y] = uvp2GP_UniG(x,y,type)
%Transfer a sepcific property(u, v, or p) to a global point location
... in a uniform grid.
    
if type=='p';
    G_x=2*x;
    G_y=2*y;
end

if  type=='v';
    G_x=2*x;
    G_y=2*y-1;
end

if  type=='u';
    G_x=2*x-1;
    G_y=2*y;

end

if  type=='None';
    G_x=2*x-1;
    G_y=2*y-1;
end

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



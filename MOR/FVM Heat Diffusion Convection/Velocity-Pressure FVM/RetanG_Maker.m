function [Grid] = RetanG_Maker(grid_Num_x,grid_Num_y,Lx,Ly,rho,k,u_T,u_B,v_L,v_R)
%% retangle grid maker
% retangle grid maker for finite volumn method.     
% Grid is a structure data including all informations about the domain
% including density, visocosity and ect.
%
% Modifications:
% 12-May-2014, WeiX, first edition 

%%
% Basic filed paramether 
Grid.rho=rho;                   % flow density
Grid.k=k;                       % Viscosity     

% Boundary condition
Grid.u_T=u_T;  Grid.v_T=0;            % Vertical & horizontal velocity of Top    lid                 
Grid.u_B=u_B;  Grid.v_B=0;            % Vertical & horizontal velocity of Bottom lid
Grid.u_L=0;  Grid.v_L=v_L;            % Vertical & horizontal velocity of Left   lid   
Grid.u_R=0;  Grid.v_R=v_R;            % Vertical & horizontal velocity of Right  lid  

% Grid information
Grid.dx=Lx/(grid_Num_x-1);
Grid.dy=Ly/(grid_Num_y-1);

Grid.Num_x_p=grid_Num_x-1; 
Grid.Num_y_p=grid_Num_y-1; 
Grid.Num_x_u=grid_Num_x;
Grid.Num_y_u=grid_Num_y-1;
Grid.Num_x_v=grid_Num_x-1;
Grid.Num_y_v=grid_Num_y;
% Grid.Num_x_GG=2*grid_Num_x-1;   % number of global grid point on x 
% Grid.Num_y_GG=2*grid_Num_y-1;
 


% uvp field 
Grid.u=zeros(Grid.Num_x_u,Grid.Num_y_u);
Grid.v=zeros(Grid.Num_x_v,Grid.Num_y_v);
Grid.p=zeros(Grid.Num_x_p,Grid.Num_y_p);

% correction coeficient for uv field
Grid.a_P_u=zeros(Grid.Num_x_u,Grid.Num_y_u);
Grid.a_P_v=zeros(Grid.Num_x_v,Grid.Num_y_v);



% % Grid.ponit=zeros(2*grid_Num_y-1,2*grid_Num_x-1);
% Grid.globalu=zeros(2*grid_Num_y-1,2*grid_Num_x-1);
% Grid.globalv=zeros(2*grid_Num_y-1,2*grid_Num_x-1);
% Grid.globalp=zeros(2*grid_Num_y-1,2*grid_Num_x-1);
% Grid.Info.dx=dx;
% Grid.Info.dy=dy;
% Grid.Info.Num_x_full=2*grid_Num_x-1;
% Grid.Info.Num_y_full=2*grid_Num_y-1;
% Grid.Info.Num_x_p=Num_x_p;
% Grid.Info.Num_y_p=Num_y_p;
% Grid.Info.Num_x_u=Num_x_u;
% Grid.Info.Num_y_u=Num_y_u;
% Grid.Info.Num_x_v=Num_x_v;
% Grid.Info.Num_y_v=Num_y_v;
end
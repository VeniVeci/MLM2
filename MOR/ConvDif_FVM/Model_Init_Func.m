function [ model ] = Model_Init_Func( Paras )
% Initialize a model with grid structure parameters and etc.
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
% Modifications:
% 6-April-2016, WeiX, first edition 

% close all

%  -----------------------Problem Parameters---------------------------
Lx=  Paras.Lx;                  % length of x axis is 1m
Ly=  Paras.Ly;                  % length of y axis is 1m
k=   Paras.k;                   % A sensitive parameter effects the stablility of solution(highter rho require more nodes for a grid)           % Thermal conductivity is 1000W/m/K     % A sensitive parameter effects the stablility of solution(highter rho require more nodes for a grid)
rho= Paras.rho;                 % flow density       % A sensitive parameter effects the stablility of solution(highter rho require more nodes for a grid)

% -----------------------Solver Parameters---------------------------
Num_node_x=Paras.Num_node_x;        % number of volumes on x direction
Num_node_y=Paras.Num_node_y;        % number of volumes on y direction
t_end=Paras.t_end;                  % total simulation time 
t_n=Paras.t_n;                      

Paras.dt=Paras.t_end/Paras.t_n;      % time step 

Paras.dx=Lx/Num_node_x;
Paras.dy=Ly/Num_node_y;
Paras.Num_node_total=Num_node_x*Num_node_y;


Node_P_index=0;

for i = 1:Num_node_x
    for j = 1:Num_node_y

        Node_P_index=Node_P_index+1;            % point index 
        
        Vol(index).CorX=dx*(1/2+(i-1));         % Coordinate X
        Vol(index).CorY=dy*(1/2+(i-1));         % Coordinate X

        Vol(index).VolW_index=Node_P_index-Num_node_y; % identify west  neighbour (in a retangular grid)
        Vol(index).VolE_index=Node_P_index+Num_node_y;
        Vol(index).VolS_index=Node_P_index-1;            % identify south neighbour (in a retangular grid)
        Vol(index).VolN_index=Node_P_index+1;    
                       
        if i==1                     % point on the left (next to boundary)  
            Vol(index).VolW_index=0;         % west neighbour is void 
        end          

        if i==Num_node_x            % point on the right
            Vol(index).VolE_index=0;      
        end         

        if j==1                     % point on the bottom
            Vol(index).VolS_index=0;      
        end         

        if j==Num_node_y            % point on the top
            Vol(index).VolN_index=0;      
        end          
        
    end
end

model.Paras=Paras;
model.Vol=Vol;


end


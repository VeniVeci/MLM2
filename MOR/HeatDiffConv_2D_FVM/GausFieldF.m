function [T,X,Y] = GausFieldF(Paras)
% Function that define velocitu in a given coordinate
%
% A prior velocity filed is governed by a given equation
%
% Modifications:
% 4-April-2016, WeiX, first edition 

% close all

dx=Paras.Lx/Paras.Num_node_x;
dy=Paras.Ly/Paras.Num_node_y;

% Num_node_total=Num_node_x*Num_node_y;


[X,Y] = meshgrid(0+dx/2:dx:Paras.Lx-dx/2,0+dy/2:dy:Paras.Lx-dy/2);

%1
Mu = [0.2 0.8];
Sigma = [0.01,0; 0,0.01];
T1 = mvnpdf([X(:) Y(:)],Mu,Sigma);
% T = reshape(T,length(x2),length(x1));
T1 = reshape(T1,Paras.Num_node_x,Paras.Num_node_y);

%2
Mu = [0.8 0.8];
Sigma = [0.01,0; 0,0.01];
T2 = 2*mvnpdf([X(:) Y(:)],Mu,Sigma);
% T = reshape(T,length(x2),length(x1));
T2 = reshape(T2,Paras.Num_node_x,Paras.Num_node_y);


%%%

%3
Mu = [0.2 0.2];
Sigma = [0.01,0; 0,0.01];
T3 = 3*mvnpdf([X(:) Y(:)],Mu,Sigma);
% T = reshape(T,length(x2),length(x1));
T3 = reshape(T3,Paras.Num_node_x,Paras.Num_node_y);

%4

Mu = [0.8 0.2];
Sigma = [0.01,0; 0,0.01];
T4 = 0*mvnpdf([X(:) Y(:)],Mu,Sigma);
% T = reshape(T,length(x2),length(x1));
T4 = reshape(T4,Paras.Num_node_x,Paras.Num_node_y);

%---------------
T=T1+T2+T3+T4;

end


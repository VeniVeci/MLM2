% kLppX_Demo
% Demonstration of kLppX 
%
% Modifications:
% WeiX, 4-11-2015, First Edition


clear
%% Parameters 
Num_data=200;
new_dim=10;

% Isomap options
options.dim_new=new_dim;                      % New dimension
options.neighborType='k';               % Type of neighbor.Choice:1)'k';Choice:2)'epsilon'
options.neighborPara=10;                % parameter for choosing the neighbor. number of neighbor for "k" type and radius for 'epsilon'
options.metric='euclidean';             % Method of measurement. Metric
options.t = 10000;

% options.kernel = 'gaussian';         % Default kernel function
options.kernel = 'linear';           % Default kernel function
% options.kernel = 'poly';             % Default kernel function
% options.kernel = 'sigmoid'; 


options.kpara = [2;0]; 
options.FullRec = 1;

%% Data Generation Swiss Roll Dataset
t1=rand(Num_data,1)*4*pi;   % Theater
t2=rand(Num_data,1)*20;
t1=sort(t1);                
X(:,1)=t1.*cos(t1);         % X
X(:,2)=t2;                  % Y
X(:,3)= t1.*sin(t1);        % Z
color = t1;                     % Color according to Theater
size = ones(Num_data,1)*10;    % Size. Constant


% Real dataset
[~,X,~,~]=Dataset_Get(100,100,5);
    
    
%% Isomaps on X    
[Z,model] = kLppX(X,options);
% [model2] = Kpca_orig(X',options);
% Z2=model2.Z';


%% Ploting
figure(1)
scatter3(X(:,1),X(:,2),X(:,3),size,color);
title('Original dataset')

figure(2)
scatter(Z(:,1),Z(:,2),size,color);
title(sprintf('Projection with %d principal components',new_dim))



% figure(3)
% scatter3(Z(:,1),Z(:,2),Z(:,3),size,color);
% title('Original dataset')



% %% KLPP by Deng Cai
% 
% % X = rand(50,10);
% % gnd = [ones(3,1);ones(15,1)*2;ones(3,1)*3;ones(15,1)*4];
% options = [];
% options.Metric = 'Euclidean';
% % options.NeighborMode = 'Supervised';
% options.NeighborMode = 'KNN';
% % options.gnd = gnd;
% options.bLDA = 1;
% W = constructW(X,options);      
% 
% options.KernelType = 'Gaussian';
% options.t = 1;
% options.Regu = 1;
% options.ReguAlpha = 0.01;
% [eigvector, eigvalue] = KLPP(W, options, X);
% 
% % feaTest = rand(5,10);
% % Ktest = constructKernel(feaTest,X,options)
% K = constructKernel(X,[],options);
% Z2 = K*eigvector;
%  
% 
% figure(3)
% scatter(Z2(:,1),Z2(:,2),size,color);
% title(sprintf('Projection with %d principal components',new_dim))
% figure(4)
% scatter3(Z2(:,1),Z2(:,2),Z2(:,3),size,color);
% title('Original dataset')
% title(sprintf('3d projecction',new_dim))



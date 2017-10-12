% LPP_Demo
% Demonstration of LPP
%
% Modifications:
% WeiX, 9-11-2015, First Edition


clear
%% Parameters 
Num_data=200;
new_dim=1;

% Isomap options
options.dim_new=new_dim;                      % New dimension
options.neighborType='k';               % Type of neighbor.Choice:1)'k';Choice:2)'epsilon'
options.neighborPara=10;                % parameter for choosing the neighbor. number of neighbor for "k" type and radius for 'epsilon'
options.metric='euclidean';             % Method of measurement. Metric
options.t = 2;                          % Parameter for heat function

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


% [~,X,~,~]=Dataset_Get(100,100,5);   %Real dataset

    
%% Isomaps on X    
% X=X';

% X = rand(50,70);
options = [];
options.Metric = 'Euclidean';
options.NeighborMode = 'KNN';
options.k = 10;
options.WeightMode = 'HeatKernel';
options.t = 1;
W = constructW(X,options);
options.PCARatio = 0.99
[eigvector, eigvalue] = LPP(W, options, X);
Z = X*eigvector;

% !!Reconstruction is not true uisng new_dim.!!


% [Z,model] = LppX(X,options);
% [model2] = Kpca_orig(X',options);
% Z2=model2.Z';

% [Z2, mapping] = lpp2(X, 2, 10, 1);

eigvector'*eigvector

%% Ploting
figure(1)
scatter3(X(:,1),X(:,2),X(:,3),size,color);
title('Original dataset')

figure(2)
scatter(Z(:,1),Z(:,2),size,color);
title(sprintf('Projection with %d principal components',new_dim))

% figure(3)
% scatter(Z2(:,1),Z2(:,2),size,color);
% title(sprintf('Projection with %d principal components',new_dim))


% figure(4)
% scatter3(Z(:,1),Z(:,2),Z(:,3),size,color);
% title('Original dataset')


%% Reconstruction
figure(4)
X_reconstruct=(eigvector'\Z')';
scatter3(X_reconstruct(:,1),X_reconstruct(:,2),X_reconstruct(:,3),size,color);
title('Reconstructed dataset')

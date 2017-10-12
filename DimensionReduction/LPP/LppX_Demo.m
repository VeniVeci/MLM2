% LppX_Demo
% Demonstration of LppX +lpp2
%
% Modifications:
% WeiX, 4-11-2015, First Edition


clear
%% Parameters 
Num_data=200;
new_dim=2;

% LPP options
options.dim_new=new_dim;                      % New dimension
options.neighborType='k';               % Type of neighbor.Choice:1)'k';Choice:2)'epsilon'
options.neighborPara=10;                % parameter for choosing the neighbor. number of neighbor for "k" type and radius for 'epsilon'
options.metric='euclidean';             % Method of measurement. Metric

options.tAuto = 1;
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

% Real dataset
% [~,X,~,~]=Dataset_Get(20,100,8);
    
%% Isomaps on X    
[Z,model] = LppX(X,options);
% [model2] = Kpca_orig(X',options);
% Z2=model2.Z';

%Compare with other codes
% [Z2, mapping] = lpp2(X, 2, 10, 1);

model.V'*model.V

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
X_reconstruct=(model.V'\Z')';
scatter3(X_reconstruct(:,1),X_reconstruct(:,2),X_reconstruct(:,3),size,color);
title('Reconstructed dataset')

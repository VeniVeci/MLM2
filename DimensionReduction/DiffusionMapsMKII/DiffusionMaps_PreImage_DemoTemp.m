% DiffusionMaps_PreImage_Demo2 MK_II
% Demonstration of PreImage solver on DiffusionMaps 
%
% Instructions:
% 1) Change the parameter section to see performance of DiffusionMaps_PreImage.
% Tips: 
% 1) options.kpara inferences the result a lot. So make sure use the proper
% one. The DiffusionMaps_AutoK function dont work quite well. Usually it
% requires for larger value for quicker desent eigenvalue, just like KPCA.
% 2) Change the input data set could lead to different result.
% 
% Modifications:
% WeiX, lost date,    first edition 
% WeiX, 11-11-2015, Minor Update
% WeiX, 21-11-2015, Minor Update

%% Initialization
clear


%% parameters
Num_data=200;
dim_new=3;

% DiffusionMap options
options.metric ='euclidean';
options.kernel ='gaussian'; 
options.kpara = 100000;             
options.kAuto=0;

options.dim_new = dim_new;              
options.t = 1;                     
options.FullRec = 0;      

options.Ztype=0;


%Diffusion PreImage options
% preoptions.type='Dw';  %'LSE','Dw' OR 'Exp'
% preoptions.para=2;
% preoptions.neighbor=10;

preoptions.type='Exp';
preoptions.neighbor=10;

% preoptions.type='LpcaI';
% preoptions.dim_new=3; % Use to stable the result but sacrefy accuracy

%-------------------------

%% Data Generation Swiss Roll Dataset
% Num_data=Num_train+Num_test;
t1=rand(Num_data,1)*4*pi;   % Theater
t2=rand(Num_data,1)*3;

t1=sort(t1);               
X(:,1)=t1.*cos(t1)*1;         % X
X(:,2)=t2*3;                  % Y
X(:,3)= t1.*sin(t1)*1;        % Z
color = t1;                     % Color according to Theater
size = ones(Num_data,1)*10;    % Size. Constant

% % Take out train & test dataset
% test_index=round(rand(Num_test,1)*Num_data);
% test_index = randperm(Num_data,Num_test);
% X_starOrigi=X(test_index,:);
% X(test_index,:)=[];


%% Main

% options.kpara=DiffusionMaps_AutoK(X,options); % Not very useful
% options.kpara=1000000;   
 
[Z,model] = DiffusionMaps(X,options);
X_star = DiffusionMaps_PreImage(Z,model,preoptions);



%% Ploting
% figure(1)
% scatter3(X(:,1),X(:,2),X(:,3),size,color);
% title('X-starOrigi dataset')
% 
% figure(2)
% scatter3(X_star(:,1),X_star(:,2),X_star(:,3),size,color);
% title('X-star dataset')
% 
% figure(3)
% plot([X_star,X]);
% title('Prediction Vs Real Curves')
% 
% figure(4)
% plot(X_star,'--')
% hold on 
% plot(X,'-')
% hold off
% title('Prediction Vs Real Curves')

figure(1)
scatter3(X(:,1),X(:,2),X(:,3),size,color);
title('Original dataset')

figure(2)
scatter(Z(:,1),Z(:,2),size,color);
title(sprintf('Projection with %d principal components',dim_new))

figure(3)
scatter3(X_star(:,1),X_star(:,2),X_star(:,3),size,color);
title(sprintf('Reconstruction of original dataset with %d principal components',dim_new))

figure(4)
plot(model.eigenvalues(1:20),'k'); %5 migh be change according to different data.
hold on 
plot(model.eigenvalues(1:dim_new),'r');
hold off
set(gca,'yscale','log');
title(sprintf('Energy Decay'))

figure(5)
plot(X_star,'--')
hold on 
plot(X,'-')
hold off
title('Prediction Vs Real Curves')





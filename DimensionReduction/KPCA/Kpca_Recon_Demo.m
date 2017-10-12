% Kpca_Recon_Demo
% Demonstration of Reconstruction of dataset X using Kpca
%
% Instructions:
% 1) Change the parameter section to see performance of Kpca_Recon_Demo.
% Tips: Larger parameter of Gaussian kernel promise a more meaning
% embedding and thus a better reconstruction. Because Larger parameter
% promises a quicker decent of eigenvalues;See model.eigenvalues
%
% Modifications:
% WeiX, Dec-3nd-2014, first edition 

%% Initialization
clear

%% parameters
Num_data=200;
dim_new=2;
              
KPCAoptions.ker='gaussian';
KPCAoptions.arg=10000;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
KPCAoptions.new_dim=dim_new;
KPCAoptions.FullRec=0;

PREoptions.type='Dw';
PREoptions.para=10;
% options.neighbor=5;

%% Data Generation Swiss Roll Dataset
% Num_data=Num_train+Num_test;
t1=rand(Num_data,1)*4*pi;   % Theater
t2=rand(Num_data,1)*20;
t1=sort(t1);                
X(:,1)=t1.*cos(t1);         % X
X(:,2)=t2;                  % Y
X(:,3)= t1.*sin(t1);        % Z
color = t1;                     % Color according to Theater
size = ones(Num_data,1)*10;    % Size. Constant

% % Take out train & test dataset
% test_index=round(rand(Num_test,1)*Num_data);
% test_index = randperm(Num_data,Num_test);
% X_starOrigi=X(test_index,:);
% X(test_index,:)=[];

%% Main
[X_star]=Kpca_Recon(X,KPCAoptions,PREoptions);


%% Ploting
figure(1)
scatter3(X(:,1),X(:,2),X(:,3),size,color);
title('Original dataset')

% figure(2)
% scatter(Z(:,1),Z(:,2),size,color);
% title(sprintf('Projection with %d principal components',dim_new))

figure(3)
scatter3(X_star(:,1),X_star(:,2),X_star(:,3),size,color);
title(sprintf('Reconstruction of original dataset with %d principal components',dim_new))

% figure(4)
% plot(real(model.eigenvalues(1:5)),'k'); %5 migh be change according to different data.
% hold on 
% plot(model.eigenvalues(1:dim_new),'r');
% hold off
% title(sprintf('Energy Decay'))

figure(5)
plot([X_star,X]);
title('Prediction Vs Real Curves')

figure(6)
plot(X_star,'--')
hold on 
plot(X,'-')
hold off
title('Prediction Vs Real Curves')
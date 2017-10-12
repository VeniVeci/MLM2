% Kpca_PreImage_Demo
% Demonstration of PreImage solver on Kpca 
%
% Instructions:
% 1) Change the parameter section to see performance of Kpca_PreImage.
% Tips: Larger parameter of Gaussian kernel promise a more meaning
% embedding and thus a better reconstruction. Because Larger parameter
% promises a quicker decent of eigenvalues;See model.eigenvalues

% Modifications:
% WeiX, lost date,    first edition 
% WeiX, Dec-2nd-2014, Add comment

%% Initialization
clear

%% parameters

% Num_train=200;
% Num_test=20;
Num_data=200;
dim_new=2;
              
options.ker='gaussian';
options.arg=1000;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
options.new_dim=dim_new;
options.FullRec=0;
options.kAuto=0;    %Use kpara=10000 give better reult than using kAuto!


% preoptions.type='Dw';
% preoptions.para=2;
% options.neighbor=10;

preoptions.type='Exp';
preoptions.neighbor=5;
preoptions.type='LpcaI';
preoptions.dim_new=3; % Use to stable the result but sacrefy accuracy


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
[Z,model] = Kpca(X,options);
X_star = Kpca_PreImage(Z,model,preoptions);


% result compare
% SSE1=sum((X_star-X).^2,2);

%% Ploting
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
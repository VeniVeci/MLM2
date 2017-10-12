% Lpca_Demo
% Demonstration of LPCA
%
% Modifications:
% WeiX, Nov-27th-2014, first edition 

clear
%% Parameters 
Num_data=1000;
dim_new=2;



%% Data Generation Swiss Roll Dataset

t1=rand(Num_data,1)*4*pi;   % Theater
t2=rand(Num_data,1)*20;
t1=sort(t1);                
X(:,1)=t1.*cos(t1);         % X
X(:,2)=t2;                  % Y
X(:,3)= t1.*sin(t1);        % Z
color = t1;                     % Color according to Theater
size = ones(Num_data,1)*10;    % Size. Constant

    
%% PCA on X    
[Z,model] = Lpca(X,dim_new);


%% Ploting
figure(1)
scatter3(X(:,1),X(:,2),X(:,3),size,color);
title('Original dataset')

figure(2)
scatter(Z(:,1),Z(:,2),size,color);
title(sprintf('Projection with %d principal components',dim_new))
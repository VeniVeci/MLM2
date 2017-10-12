% Lpca_PreImage_Demo
% Demonstration of Lpca_PreImage
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
[X_star]=Lpca_PreImage(Z,model);


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
plot(model.eigenvalues,'k');
hold on 
plot(model.eigenvalues(1:dim_new),'r');
hold off
title(sprintf('Energy Decay'))

figure(5)
plot([X_star,X]);
title('Prediction Vs Real Curves')

figure(6)
plot(X_star,'--')
hold on 
plot(X,'-')
hold off
title('Prediction Vs Real Curves')


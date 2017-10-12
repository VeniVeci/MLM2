% Kpca_Demo
% Demonstration of KPCA
%
% Instructions:
% 1) Change the parameter section to see performance of Kpca_PreImage.
% Tips: Larger parameter of Gaussian kernel promise a more meaning
% embedding and thus a better representation. Because Larger parameter
% promises a quicker decent of eigenvalues;See model.eigenvalues
%
% Modifications:
% WeiX, Dec-1th-2014, first edition 

clear
%% Parameters 
Num_data=1000;
new_dim=2;

options.ker='gaussian';
options.arg=100;         %using 10 to show a fail case
options.new_dim=new_dim;
options.FullRec=0;
options.kAuto=0;

%% Data Generation Swiss Roll Dataset

t1=rand(Num_data,1)*4*pi;   % Theater
t2=rand(Num_data,1)*10;
t1=sort(t1);                
X(:,1)=t1.*cos(t1);         % X
X(:,2)=t2;                  % Y
X(:,3)= t1.*sin(t1);        % Z
color = t1;                     % Color according to Theater
size = ones(Num_data,1)*10;    % Size. Constant

    
%% KPCA on X    
[Z,model] = Kpca(X,options);
% [model2] = Kpca_orig(X',options);
% Z2=model2.Z';


%% Ploting
figure(1)
scatter3(X(:,1),X(:,2),X(:,3),size,color);
title('Original dataset')

figure(2)
scatter(Z(:,1),Z(:,2),size,color);
title(sprintf('Projection with %d principal components',new_dim))





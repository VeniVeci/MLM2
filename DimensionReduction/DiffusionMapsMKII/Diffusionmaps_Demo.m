% Diffusionmaps_Demo
% Demonstration of Diffusionmaps
%
% Modifications:
% WeiX, Dec-5th-2014, first edition 

clear
%% Parameters 
Num_data=1000;

% DiffusionMap options
options.metric ='euclidean';
options.kernel ='gaussian'; 
options.kpara = 10000;    
options.kAuto=1;

options.dim_new = 5;              
options.t = 1;                     
options.FullRec = 0;      


%% Data Generation Swiss Roll Dataset
t1=rand(Num_data,1)*4*pi;   % Theater
t2=rand(Num_data,1)*20;     % change 20 to 10 will show different
t1=sort(t1);                
X(:,1)=t1.*cos(t1);         % X
X(:,2)=t2;                  % Y
X(:,3)= t1.*sin(t1);        % Z
color = t1;                     % Color according to Theater
size = ones(Num_data,1)*10;    % Size. Constant

    
%% Diffusionmaps on X    
% options.kpara=DiffusionMaps_AutoK(X,options); % Not very useful
% options.kpara=100;   
 
[Z,model] = DiffusionMaps(X,options);

%% Ploting
figure(1)
scatter3(X(:,1),X(:,2),X(:,3),size,color);
title('Original dataset')

figure(2)
scatter(Z(:,2),Z(:,3),size,color);
title(sprintf('Projection with %d principal components',options.dim_new))


%%

% for i=1:options.dim_new
%     p1_recover=Z(1,1:i)*model.L(:,1:i)';
%     p_recover(i,:)=p1_recover;        
% end
% 
% pr_sum=sum(p_recover);
% figure
% plot(pr_sum);
% 
% p_orig=model.P;
% 
% 
% p_recover=p_recover';   
% figure
% plot(p_recover);
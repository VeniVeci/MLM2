% Dist2pos_Demo
% Demonstration of distances to position(coordinate) conversion.
%
% Modifications:
% WeiX, lost date,    first edition 
% WeiX, Dec-2nd-2014, Add comment

%% Initialization
clear

%% parameters

Num_train=40;
Num_test=200;
Dim=39;                 % When Dim >=Num_train, LSE algorithm would fail. Adjust this parameter to see more. For Random 'Data Generation Random' data below

options.type='Dw';
options.para=2;
options.neighbor=10;


%% Data Generation Random data
X=rand(Num_train,Dim)*10;
X_starOrigi=rand(Num_test,Dim)*10;

% %% Data Generation Swiss Roll Dataset
% Num_data=Num_train+Num_test;
% t1=rand(Num_data,1)*4*pi;   % Theater
% t2=rand(Num_data,1)*20;
% t1=sort(t1);                
% X(:,1)=t1.*cos(t1);         % X
% X(:,2)=t2;                  % Y
% X(:,3)= t1.*sin(t1);        % Z
% color = t1;                     % Color according to Theater
% size = ones(Num_data,1)*10;    % Size. Constant
% 
% % Take out train & test dataset
% % test_index=round(rand(Num_test,1)*Num_data);
% test_index = randperm(Num_data,Num_test);
% X_starOrigi=X(test_index,:);
% X(test_index,:)=[];



%% Main
for i = 1:Num_test
   %calculate distance
   disti=pdist2(X,X_starOrigi(i,:));
   
   %Recover using distances
   X_star1(i,:)= Dist2pos(X,disti,options);
   
end

% result compare
SSE1=sum((X_star1-X_starOrigi).^2,2);


%% Ploting
figure(1)
scatter3(X_starOrigi(:,1),X_starOrigi(:,2),X_starOrigi(:,3));
title('X-starOrigi dataset')

figure(2)
scatter3(X_star1(:,1),X_star1(:,2),X_star1(:,3));
title('X-star dataset')

figure(3)
plot(SSE1)
title('Square sum error of each test point')



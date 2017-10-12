% Kpca_ReconAnal
% Demonstration of Reconstruction & Analysis of dataset X using Kpca
%
% Modifications:
% WeiX, Dec-3nd-2014, first edition 

%% Initialization
clear

%% parameters
Num_data=400;
              
KPCAoptions.ker='gaussian';
KPCAoptions.arg=10000;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
KPCAoptions.new_dim=2;
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
[X_star,Anal,model]=Kpca_ReconAnal(X,KPCAoptions,PREoptions);


%% Ploting
figure(1)
scatter3(X(:,1),X(:,2),X(:,3),size,color);
title('X-starOrigi dataset')

figure(2)
scatter3(X_star(:,1),X_star(:,2),X_star(:,3),size,color);
title('X-star dataset')

figure(3)
plot(X_star(:,1),'r--');
hold on 
plot(X(:,1),'k');
hold off
title('Dimension 1st')

figure(4)
plot(X_star(:,2),'r--');
hold on 
plot(X(:,2),'k');
hold off
title('Dimension 2nd')

figure(5)
plot(X_star(:,3),'r--');
hold on 
plot(X(:,3),'k');
hold off
title('Dimension 3rd')

Fix_RateAbsDiff=Anal.Distance_RateAbsDiff.*(Anal.Distance_RateAbsDiff<1)+(Anal.Distance_RateAbsDiff>1);
figure(6)
mesh(Fix_RateAbsDiff);
zlim([0,1])
title('Distance different rate. Exceeding 1 mean failure')

figure(7)
mesh(Anal.Distance_orig);
title('Original distance map')

figure(8)
mesh(Anal.Distance_star);
title('Reconstructed distance map')

figure(9)
plotmatrix(X,X_star)
title('Original VS Reconstructed')

% DiffusionMaps_PreImage_Demo MK_II
% Demonstration of PreImage solver on DiffusionMaps 
%
%
% Modifications:
% WeiX, 22-12-2015, first edition 
clear 

Num=500;
Doptions.para=0.03;
dim_new=2;


% Type='SwissRoll';
% Type='SwissHole';
% Type='CornerPlanes';
% Type='PuncturedSphere';
% Type='TwinPeaks';
% Type='3DClusters';
% Type= 'ToroidalHelix';
% Type= 'Gaussian';
Type= 'Spiral';

[Data] = Data_GeneratorT(Num,Type,Doptions);
% [Data] = Data_GeneratorM(Num,Type,Doptions);

% DiffusionMap options
options.metric ='euclidean';
options.kernel ='gaussian'; 
options.kpara =1000000;     %10000 works much better than auto kpara in Type='Gaussian'        
options.kAuto=0;

options.dim_new = dim_new;              
options.t = 1;                     
options.FullRec = 1;      

options.Ztype=0;

%Diffusion PreImage options
% preoptions.type='Dw';  %'LSE','Dw' OR 'Exp'
% preoptions.para=5;
% preoptions.neighbor=10;

preoptions.type='Exp';
preoptions.neighbor=10;
% preoptions.type='LpcaI';
% preoptions.dim_new=3; % Use to stable the result but sacrefy accuracy



X=Data.Y;
[Z,model] = DiffusionMaps(X,options);
X_star = DiffusionMaps_PreImage(Z,model,preoptions);


squErr=(X_star-X).^2;
squSumErr=sum(squErr(:))
averageSSE=sum(squErr(:))/Num
relativeErr=sum(squErr,2)./sum(X_star,2);
averageRe=sum(relativeErr)/Num


% figure(1)
% scatter3(X(:,1),X(:,2),X(:,3),size,color);
% title('Original dataset')
% 
% figure(2)
% scatter(Z(:,1),Z(:,2),size,color);
% title(sprintf('Projection with %d principal components',dim_new))
% 
% figure(3)
% scatter3(X_star(:,1),X_star(:,2),X_star(:,3),size,color);
% title(sprintf('Reconstruction of original dataset with %d principal components',dim_new))
% 
% figure(4)
% plot(model.eigenvalues(1:20),'k'); %5 migh be change according to different data.
% hold on 
% plot(model.eigenvalues(1:dim_new),'r');
% hold off
% set(gca,'yscale','log');
% title(sprintf('Energy Decay'))
% 
% figure(5)
% plot(X_star,'--')
% hold on 
% plot(X,'-')
% hold off
% title('Prediction Vs Real Curves')




figure(1)
scatter3(Data.Y(:,1),Data.Y(:,2),Data.Y(:,3),Data.SizeVector,Data.ColorVector)
title( sprintf( 'Original data %d points',Num)) 

figure(2)
scatter3(X_star(:,1),X_star(:,2),X_star(:,3),Data.SizeVector,Data.ColorVector)
title( sprintf( 'Diffusion maps reconstructin with %d dimension data',dim_new)) 

figure(3)
scatter(Z(:,1),Z(:,2),Data.SizeVector,Data.ColorVector)
title( sprintf( 'Diffusion maps %d-D representation',dim_new)) 


figure(4)
plot(X_star,'--')
hold on 
plot(X,'-')
hold off
title('Prediction Vs Real Curves')

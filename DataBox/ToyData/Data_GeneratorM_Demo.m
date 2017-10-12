%% Data_GeneratorM_Demo file
%  Data_GeneratorM Demo file
%
%  Modification
%  WeiX Nov 25th 2014 First Edition
clear 

Num=1000;
options.para=0.03;

Type='SwissRoll';
Type='SwissHole';
Type='CornerPlanes';
Type='PuncturedSphere';
Type='TwinPeaks';
Type='3DClusters';
Type= 'ToroidalHelix';
Type= 'Gaussian';
Type= 'Spiral';


[Data] = Data_GeneratorM(Num,Type,options);

figure(1)
scatter3(Data.Y(:,1),Data.Y(:,2),Data.Y(:,3),Data.SizeVector,Data.ColorVector)

% figure(2)
% scatter(Data.X(:,1),Data.X(:,2),Data.SizeVector,Data.ColorVector)

figure(2)
scatter3(Data.Y(:,1),Data.Y(:,2),Data.Y(:,3));
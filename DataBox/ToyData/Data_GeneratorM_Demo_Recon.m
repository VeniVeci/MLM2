% Data_GenerationM_Demo_Recon
% Demonstration of reconstruction using Dimension reductions
%
% Instructions:
% 1) Change the parameter section to see performance of Kpca_PreImage.
%
% Modifications:
% WeiX, Dec-8nd-2014, First Edition
%

clear 
%% Data Parameter 
Num=500;
dataoptions.para=5;
Type='SwissRoll';
% Type='SwissHole';
% Type='CornerPlanes';
% Type='PuncturedSphere';
% Type='TwinPeaks';
% Type='3DClusters';
% Type= 'ToroidalHelix';
% Type= 'Gaussian';

%% Common parameter
dim_new=2;

%% Kpca Parameter 
kpca.active=1;

kpca.options.ker='gaussian';
kpca.options.arg=10000;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
kpca.options.new_dim=dim_new;
kpca.options.FullRec=0;

% kpca PreImage options
kpca.preoptions.type='Dw';
kpca.preoptions.para=10;
% options.neighbor=5;

%% Isomap Parameter
isomaps.active=1;

isomaps.options.dim_new=dim_new;                % New dimension
isomaps.options.neighborType='k';               % Type of neighbor.Choice:1)'k';Choice:2)'epsilon'
isomaps.options.neighborPara=5;                 % parameter for choosing the neighbor. number of neighbor for "k" type and radius for 'epsilon'
isomaps.options.metric='euclidean';             % Method of measurement. Metric

% Isomap PreImage options
isomaps.preoptions.neighborType='k';    % Type of neighbor of new point. Choice:1)'k';Choice:2)'epsilon'
isomaps.preoptions.neighborPara=5;      % Parameter of neighbor of new point
isomaps.preoptions.type='Dw';           % Type of distance to coordinate method. Distance weight/Least square estimate
isomaps.preoptions.para=1;              % Parameter of distance to coordinate recover method

%% Diffusiomaps Parameter
diffusionmaps.active=1;

diffusionmaps.options.metric ='euclidean';
diffusionmaps.options.kernel ='gaussian'; 
diffusionmaps.options.kpara = 10000;             
diffusionmaps.options.dim_new = dim_new;              
diffusionmaps.options.t = 1;     % 10 for 'Gaussian' type data for interesting result                 
diffusionmaps.options.FullRec = 0;      

% Diffusion PreImage options
diffusionmaps.preoptions.type='Dw';  %'LSE' OR 'Dw'
diffusionmaps.preoptions.para=2;
diffusionmaps.preoptions.neighbor=5;

%% Data Generation
[Data] = Data_GeneratorM(Num,Type,dataoptions);
X=Data.Y;


%% Main 
if kpca.active==1
    [kpca.Z,model] = Kpca(X,kpca.options);
    kpca.X_star = Kpca_PreImage(kpca.Z,model,kpca.preoptions);
end

if isomaps.active==1
    [isomaps.Z,model] = Isomaps(X,isomaps.options);
    isomaps.X_star = Isomaps_PreImage(isomaps.Z,model,isomaps.preoptions);
end

if diffusionmaps.active==1
%     [ diffusionmaps.options.kpara ] = DiffusionMaps_AutoK(X,diffusionmaps.options);
    [diffusionmaps.Z,model] = DiffusionMaps(X,diffusionmaps.options);
    diffusionmaps.X_star = DiffusionMaps_PreImage(diffusionmaps.Z,model,diffusionmaps.preoptions);
end

%% Plot
figure(1)
scatter3(X(:,1),X(:,2),X(:,3),Data.SizeVector,Data.ColorVector);
title('Original dataset')


if kpca.active==1
    figure(2)
    scatter3(kpca.X_star(:,1),kpca.X_star(:,2),kpca.X_star(:,3),Data.SizeVector,Data.ColorVector);
    title(sprintf('Kpca Reconstruction with %d principal components',dim_new))

    figure(3)
    scatter(kpca.Z(:,1),kpca.Z(:,2),Data.SizeVector,Data.ColorVector);
    title(sprintf('Kpca Projection with %d principal components',dim_new))
end

if isomaps.active==1
    figure(4)
    scatter3(isomaps.X_star(:,1),isomaps.X_star(:,2),isomaps.X_star(:,3),Data.SizeVector,Data.ColorVector);
    title(sprintf('Isomaps Reconstruction with %d principal components',dim_new))
    figure(5)
    scatter(isomaps.Z(:,1),isomaps.Z(:,2),Data.SizeVector,Data.ColorVector);
    title(sprintf('Isomaps Projection with %d principal components',dim_new))
end

if diffusionmaps.active==1
    figure(6)
    scatter3(diffusionmaps.X_star(:,1),diffusionmaps.X_star(:,2),diffusionmaps.X_star(:,3),Data.SizeVector,Data.ColorVector);
    title(sprintf('Diffusionmaps Reconstruction with %d principal components',dim_new))

    figure(7)
    scatter(diffusionmaps.Z(:,1),diffusionmaps.Z(:,2),Data.SizeVector,Data.ColorVector);
    title(sprintf('Diffusionmaps Projection with %d principal components',dim_new))
end































figure(1)
scatter3(Data.Y(:,1),Data.Y(:,2),Data.Y(:,3),Data.SizeVector,Data.ColorVector)
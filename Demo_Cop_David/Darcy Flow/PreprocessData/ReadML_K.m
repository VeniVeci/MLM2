%% SCRIPT
% Manifold learning on the input data read
% Export decomposition on filed input data.


% currentFolder=pwd;
% DataFolder='./PreprocessData/';
% cd(DataFolder)
% % figList= dir('*.mat');
% cd(currentFolder)


clear 

load('K_kPCA.mat')
eigenvalues=model.eigenvalues.*(model.eigenvalues>0);
eigenPercentage(1,:)=cumsum(eigenvalues)/sum(eigenvalues);

load('K_Isomaps.mat')
eigenvalues=model.eigenvalues.*(model.eigenvalues>0);
eigenPercentage(2,:)=cumsum(eigenvalues)/sum(eigenvalues);

load('K_DiffusionMaps.mat')
eigenvalues=model.eigenvalues.*(model.eigenvalues>0);
eigenPercentage(3,:)=cumsum(eigenvalues)/sum(eigenvalues);

load('K_PCA.mat')
eigenvalues=model.eigenvalues.*(model.eigenvalues>0);
eigenPercentage(4,:)=cumsum(eigenvalues)/sum(eigenvalues);

% load ('./RawData/Eigenmodes/lambda_0_3/S_elliptic.mat');
load ('S_elliptic.mat');
eigenvalues=[diag(S);0];
cum4All=cumsum(eigenvalues)/sum(eigenvalues);
eigenPercentage(5,:)=cum4All(1:1000,:);



eigenPercentage=eigenPercentage';

plot(eigenPercentage,'DisplayName','eigenPercentage')
legend('KPCA','Isomaps','DiffusionMaps','PCA','PCA on Sigma')
title('Acumulated eigenvalue percentage directly on K');


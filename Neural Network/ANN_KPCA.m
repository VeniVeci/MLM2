function [Y_star_kpca,time]=ANN_KPCA(X,Y,X_star,options,kfold)
% ANN prediction with PCA via KPCA
%  
% Synopsis:
% [Y_star_kpca,time]=ANN_KPCA(X,Y,X_star,options,kfold)
%
% Description:
% Use KPCA to truncate the dataset in output dataset Y and then use the
% truncated data to feed to ANN to make predictions.
% 
% steps for the program:
% 
% Input:
%  X        [dimension_X x samples] indicates Training data X.
%  Y        [dimension_X x samples] indicates Training data Y.
%  X_star   [dimension_Y x test points] indicates conditions of prediction Y_star
%  options [struct] Decribes kernel type for KPCA 
%   .ker [string] Kernel identifier ('linear','poly','rbf','sigmoid'); 
%   .arg [1 x 2] kernel argument; 
%   .new_dim [1x1] dimensions in used (number of used principal 
%     components);
%  Kfold    [1 x 1] value unmber indicates fold used for crossvalidation
%
% Output:
%  Y_star      [dimension_Y x test points] mean values of predictions
%  Y_var_star  [dimension_Y x test points] variance of preditions
%  time            [1 x 1]   time used for the algorithm
%
% Example:
%
% See also 
%  Kernel, Kpca. 
% 
% About: 
% 
% Modifications:
% 19-jul-2013, WeiX, first edition 


[~,np_train]=size(X);
[~,np_test] =size(X_star);
[Dim_Y,~]=size(Y);

tic;
model2 = Kpca(Y,options);

Z_star_kpca= ANN_prediction_crossvalidation(X,model2.Z,X_star,kfold);

Y_star_kpca = Kpca_PreImage(Z_star_kpca,model2);

time=toc;

return

function [Y_star_SVD,time]=ANN_SVD(X,Y,X_star,new_dim,kfold)
% ANN prediction with with PCA via SVD
%  
% Synopsis:
% [Y_star_kpca,time]=ANN_KPCA(X,Y,X_star,options,kfold)
%
% Description:
% Use PCA to truncate the dataset in output dataset Y and then use the
% truncated data to feed to ANN to make predictions.
% 
% steps for the program:
% 
% Input:
%  X        [dimension_X x samples] indicates Training data X.
%  Y        [dimension_X x samples] indicates Training data Y.
%  X_star   [dimension_Y x test points] indicates conditions of prediction Y_star
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

%---------------------------
%change the formate of data
Y=Y';
X=X';
X_star=X_star';
%---------------------------
%Testing

tic;
[Z_SVD,Key] = DataReduc_SVD(Y,new_dim);
[Z_star_SVD]= ANN_prediction_crossvalidation(X',Z_SVD',X_star',kfold);

Z_star_SVD=Z_star_SVD';
Y_star_SVD=Z_star_SVD*Key;
Y_star_SVD=Y_star_SVD';

time=toc;

return

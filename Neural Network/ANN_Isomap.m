function [Y_star,time]=ANN_Isomap(X,Y,X_star,options_Isomap,options_Remap,kfold)
% Neural Network Regression Prediction with isomap 
%
% Synopsis:
%  [Y_star_kpca,Yvar_star_kpca,time]=GPR_Isomap(X,Y,X_star,options)
%
% Description:
% Use isomap to truncate the dataset in output dataset Y and then use the
% truncated data to feed to gaussian process regression to make
% predictions.
% 
% steps for the program:
% 
% Input:
%  X        [dimension_X x samples] indicates Training data X.
%  Y        [dimension_X x samples] indicates Training data Y.
%  X_star   [test points x dimension_Y] indicates input of prediction Y_star
%  options  [struct] parameters for isomap and its preimage solution
%   .dim_new [1x1] dimensions in used 
%   .neighbor[1x1]= number of neighbor points when calculate the distance map
%   .d2p_method indicates the method when reconstruct coordinates using
%    distances.=either "Dw"(Distance weight) or"LSE"(Least square estimate.
%   .d2p_Dwpara [1x1] indicates the power for weight when using "Dw" method
%   .d2p_points [1x1] indicates the number of closest points when
%    reconstruct using distances
% 
% Output:
% Y_star      [dimension_Y x test points] mean values of predictions
% time        [1 x 1]   time used for the algorithm
%
% Example:
%
% See also 
%  Isomap, MDS. 
% 
% About: 
% 
% Modifications:
% 27-Aug-2014, WeiX, Create
%
% options = struct('ker','rbf','arg',para_kpca,'new_dim',new_dim);

tic;

Y=Y';

[Z,model] = Isomap(Y,options_Isomap);

Z=Z';

Z_star= ANN_prediction_crossvalidation(X,Z,X_star,kfold);

Z_star=Z_star';

Y_star = Isomap_PreImage(Z_star,model,options_Remap);

Y_star=Y_star';

time=toc;

return

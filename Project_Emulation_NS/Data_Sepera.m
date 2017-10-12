function [ X_train,X_test ] = Data_Sepera(X,num_train,num_test,options)
%Data seperarion for training and testing
%  
% Synopsis:
% [ X_train,X_test ] = Data_Sepera( X )
%
% Description:
% 
% steps for the program:
% 
% Input:
%  X        [samples_x dimension_X] indicates Training input data X.
% options.order
% 
% Output:
% X_train      [Train points x dimension_X] mean values of predictions
% X_test       [Test points x  dimension_X] variance of preditions
%
% Example:
%
% See also 
%  Kernel, Kpca. 
% 
% About: 
% 
% Modifications:
% 18-Nov-2014, WeiX, first edition 

[num,dim]=size(X);

if ~isfield(options,'order'), options.order ='normal'; end    

if nargin <=1
    num_train=round(num/5);
    num_test=point-num_train;
end

if (num_train+num_test)>num
    error('Training number and testing number overlap ')
end

switch options.order
    case 'normal'
        X_train=X(1:num_train,:);
        X_test=X(end-num_test+1:end,:);
    otherwise 
        error('Order method not developed Yet!')
end

            
   

end


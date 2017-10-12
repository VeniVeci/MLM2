function [X_star]=Lpca_Recon(X,dim_new,options)
% Lpca Reconstruction via preimage solver
%
% Synopsis:
% [X_star]=Lpca_Recon(X,X_test,new_dim)
%
% Description:
% PCA is run on X to achieve reduction dataset Z. Z is then used to 
% reconstruct X using bases(principal direction) offered by PCA on X.
% # out of sample is not a concern here.
%
%
% Input:
%  X  [sample x dimension] indicates Training data.
%  dim_new [1 x 1] value number indicates the principal components that is
%  used in the dimension reduction process.
% 
%  X_star [sample x dimension] Thhe reconstructed original dataset X from Z
% 
% Example:
%
% See also 
% Lpca, Lpca_PreImage
% 
% About: 
% 
% Modifications:
% 18-Nov-2013, WeiX, first edition 
% WeiX, Nov-27th-2014, 2nd version & structure update

%% Initialization and Parameters
if nargin < 3, options = []; end
if ~isfield(options,'FullRec'), options.FullRec = 0; end                 % Default output information;

%% Reconstruction
[Z,model] = Lpca(X,dim_new,options);
X_star=Z*model.Key;


%-------------------------------------------------------------------------
% % Version II
% [Num_train,Dim_X]=size(X_train);
% [Num_test,Dim_X]=size(X_test);
% 
% [Z,Key] = DataReduc_SVD(X_train,dim_new);
% Z_star=X_test*Key';
% X_test_star=Z_star*Key;

%-------------------------------------------------------------------------
% % Version I
% for j = 1:Num_test
%     
%     % carrying out KPCA on X+X_test(i) 
%     X=[X_train;X_test(j,:)];
% 
%     [Z,Key] = DataReduc_SVD(X,new_dim);
%     Z_star=Z(end,:);
%     X_test_star(j,:)=Z_star*Key;
%     
% end
%-------------------------------------------------------------------------

        
return

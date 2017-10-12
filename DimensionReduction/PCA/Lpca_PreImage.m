function [X_star]=Lpca_PreImage(Z_star,model)
% Lpca_PreImage solver
%
% Synopsis:
% [X_star]=Lpca_PreImage(Z_star,model)
%
% Description:
% Find original coordinate of Z_star using the principal directions(Key).
% #2nd method is to use loca linear interpretation. Not yet developed.
%
% Input:
%  Z_star   [sample x dimension]  data point in reduced space.
%  model    [structure]           Information about the original data and model
%
% Output:
%  X_test   [sample x dimension]  corresponding data point in original space
%
% Example:
%
% See also 
% Lpca. 
% 
% Modifications:
% 18-Nov-2013, WeiX, first edition 
% WeiX, Nov-27th-2014, 2nd version update.

[Num_Zstar,Dim_Zstar]=size(Z_star);

X_star=Z_star*model.Key(1:Dim_Zstar,:);


%-------------------------------------------------------------------------
% % Version II
% [Num_train,Dim_X]=size(X_train);
% [Num_test,Dim_X]=size(X_test);
% 
% [Z,Key] = DataReduc_SVD(X_train,new_dim);
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

        
end

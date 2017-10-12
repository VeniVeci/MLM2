function [Z,model] = Lpca(X,dim_new,options) % data_replacement(dimension reduced)*Key ~= data 
% Data dimension reduction using SVD methods
%  
% Synopsis:
% [Z,Key] = DataReduc_SVD(data,dim_new)
%
% Description:
% This function apple PCA by doing SVD to data
% z*Key ~= data 
% 
% steps for the program:
% 
% 
% Input:
% X                         % [sample x dimension] indicates Training data.
% dim_new                   % [1 x 1] value indicates the principal components used
%                           for PCA. dim_new <= dimension !!
% 
% Output:
% Z [sample x dim_new]      projections of data via PCA.
% 
% model [structure]         % Information about the original data and model
%      .Key                 % [dim_new x dimension] Principal directions used to recover
%                             the original data. !Z*Key ~= data!
%      .DR_method           % Dimension reduction method
%      .dim_new=dim_new;    % New dimension
%      .X=X;                % Original dataset in original space
%      .Key=Key;            % Recover Key. Principal directions.
%      .energy=diag(S);     % Energy series. Eigenvalue of covariance matrix
%
% Example:
%
% See also:
% Lpca_PreImage
% 
% About: 
% 
% Modifications:
% 19-jul-2013, WeiX, first edition 
% WeiX, Nov-27th-2014, 2nd version update.

%% Initialization and Parameters
% [num,dim]=size(X);
% if nargin < 2, options = []; else options=c2s(options); end
if nargin < 3, options = []; end
if ~isfield(options,'FullRec'), options.FullRec = 0; end                 % Default output information;

%% PCA process via SVD
[U,S,V]=svd(X);                     % X=U*S*V'
% [U,S,V]=svd(X,'econ');                % for data that is large

U_temp=U(:,1:dim_new);
S_temp=S(1:dim_new,1:dim_new);      
V_temp=V(:,1:dim_new);              % X~=U_temp*S_temp*V_temp'

Z=U_temp*S_temp;                    % Projections
Key=V_temp';                        % Principal directions
% Z*Key ~= data 
% new data (predited one) can also projected to the original formate by
% times Key (see instruction inference)

%% Recording & Output
model.DR_method='LPCA';
model.dim_new=dim_new;
model.X=X;
model.Key=Key;
model.eigenvalues=diag(S); % also known as erengy

% Full Information. Conditional output. (For Further Research without recalculation. Mass Memoey required)
if options.FullRec == 1
    model.U=U;              % Full U Eigenvectors 
    model.S=S;              % Full Lambda Eigenvalues
    model.V=V;              % Full V Eigenvectors 
end

end

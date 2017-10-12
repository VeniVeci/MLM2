function [X_star,model] = LpcaI(X,Z,Z_star,options) 
% Local (linear) PCA interpretation
% Use neighbor to build PCA basis and LSE to interpret corresponding
% coordinate to build mapping. Designed to approximate inverse mapping for
% nonlinear dimension reduction technique.
%  
% Synopsis:
% [X_star,model] = LpcaI(X,Z,Z_star,options) 
%
% Input:
% X                         % [sample x dimension] observation in original space
% Z                         % [sample x dim_fea] corresponding feature observations
% Z_star                    % [num x dim_fea] feature observations to be
% options                   % parameters
%
% Output:
% X_star                    % [num x dimension] recovered observation in original space 
% 
% model [structure]         % Information about the original data and model
%      .options             % Function parameters
%
% Example:
%
% See also:
% Lpca_PreImage & Lpca
% 
% About: 
% 
% Modifications:
% WeiX, 1-11-2016, DayI Create 

%% Initialization and Parameters
[num_Z_star,DimF]=size(Z_star);     %Dimension of feature space
% if nargin < 2, options = []; else options=c2s(options); end
if nargin < 4, options = []; end
if ~isfield(options,'dim'), options.dim_new = DimF; end                     % Default output information;
if ~isfield(options,'InMethod'), options.InMethod = 'LSE'; end                     % Default Interpretation method;
if ~isfield(options,'FullRec'), options.FullRec = 0; end                    % Default output information;


%% PCA process via SVD on X
dim_new=options.dim_new;
[U,S,V]=svd(X);                     % X=U*S*V'

U_temp=U(:,1:dim_new);
S_temp=S(1:dim_new,1:dim_new);      
V_temp=V(:,1:dim_new);              % X~=U_temp*S_temp*V_temp'

Y=U_temp*S_temp;                    % Projections
% Key=V_temp';                        % Principal directions
% Z*Key ~= data 
% new data (predited one) can also projected to the original formate by
% times Key (see instruction inference)

%% Interprete new coordinate
switch options.InMethod
    case 'LSE'
        A=Y\Z;              %Linear approximation Z to Y via LSE
        Y_star=A*Z_star;
        X_star=Y_star*V_temp';
    otherwise 
        error('No such Interpret Method')
end
    

%% Recording & Output
model.options=options;

% Full Information. Conditional output. (For Further Research without recalculation. Mass Memoey required)
if options.FullRec == 1
    model.U=U;              % Full U Eigenvectors 
    model.S=S;              % Full Lambda Eigenvalues
    model.V=V;              % Full V Eigenvectors 
end



end

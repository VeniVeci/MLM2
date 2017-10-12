function [Z,model] = DiffusionMaps(X,options)
% Description:
% Diffusion map Dimension Reduction Function
%
% Synopsis:
% [Z,model] = Diffusionmap(X)
% [Z,model] = Diffusionmap(X,option)
%
% Input: 
%       X              % Original dataset [samples X dimension ] matrix
%       option         % Options of the diffusion map process
%             .kernel  % Type of kernel function
%             .kpara   % Parameter of the kernel function
%             .metric  % Method of measurement. Metric
%             .t       % optional time parameter in diffusion map 
%                          (default: model with multiscale geometry)    
%             .FullRec % Record all output information flag. 0 for false
%                          and 1 for ture.
%
% Output: 
%       Z              % new coordinates system [new dimension x samples] 
%       model          % preserving the Diffusion map process information and
%                      parameters
%
% About: 
%  Modification
%  WeiX, 3-11-2014, First Edition

%
%
%% Initialization and Parameters
% [num,dim]=size(X);
if nargin < 2, options = []; end
if ~isfield(options,'metric'), options.metric ='euclidean'; end          % Default metric to measure the distance;
if ~isfield(options,'kernel'), options.kernel = 'gaussian'; end          % Default kernel function
if ~isfield(options,'kpara'),  options.kpara = 1000; end                 % Default kernel parameter
if ~isfield(options,'dim_new'),options.dim_new = 2; end                  % Default new dimension=3;
if ~isfield(options,'t'), options.t = 1; end                             % Default Diffusion times;
if ~isfield(options,'FullRec'), options.FullRec = 0; end                 % Default output information;

dim_new=options.dim_new;

%% Main
% Calculating distances
Distance =pdist2(X,X,options.metric);

% Calculating the Kernel Matrix
switch options.kernel
    case 'gaussian'
        K = exp(-Distance.^2/options.kpara);            
    otherwise
        error('Error: Undefined type of kernel function.');    
end

d=sum(K,2);                             % row sum
D=diag(d);                              % Degree Matrix. It is a diagonal matrix. 
L=D^(-1)*K;                             % The Laplacian Matrix
Ln=D^(-0.5)*K*D^(-0.5);                 % The Normalized Laplacian Matrix

%-----------------------------------------
% Eigen Decomposition of The Laplacian Matrix
[V,Lamda] = eig(L);  
[Lamda,I] = sort(diag(Lamda),'descend'); % sort eigenvalues in descending order
V = V(:, I);                             % sort eigenvector in corresponding order to its eigenvalue
Lamda=diag(Lamda);

%-------------------------------------------------
% Detection. 
% This part is made to ensure DiffusionMap is carried out safely
% under the limit of computer. Situation happens when K(i,j)=0.
if abs(mean(V(:,1))-V(1,1))>1e-3  % then each element in V are assumed to be the same 
    error('parameter "options.kpara" in Diffusion Map is not suitable(Probably too small). Please Change.') 
end

if ~isreal(V)
   warning('Eigenvector V of Laplacian Matrix L is not real number. Result might not be accurate. Parameter "options.kpara" is adviced to change. This could be normal as the few last V tend to be complex number vector') 
end

% if sum(K(:,1))>length(K(:,1))/3
%     warning('Kernel Matrix varies very limit. Parameter "options.kpara" in Diffusion Map is not suitable(Usually too large). Please Change. ') 
% end

if (max(K(:,1))-min(K(:,1)))<0.9    
    warning('Kernel Matrix varies very limit. Parameter "options.kpara" in Diffusion Map is not suitable(Usually too large). Please Change. ') 
end
    
%------------------------------------------
% Eigen Decomposition of The Normalized Laplacian Matrix
[Vn,Lamdan] = eig(Ln);  
[Lamdan,I] = sort(diag(Lamdan),'descend');
Vn = Vn(:, I);
Lamdan=diag(Lamdan);                     % Lamdan must = Lamda

% ------------------------------------------
% % Truncation/Dimension Reduction
% V=V(:,1:Dim_new);
% Lamda=Lamda(1:Dim_new,1:Dim_new);
% 
% Vn=Vn(:,1:Dim_new);
% Lamdan=Lamdan(1:Dim_new,1:Dim_new);

% -----------------------------------------
% Relation Matrix between V and Vn. It is a diagonal matrix
C=D^(-0.5)*Vn./V;
C=diag(C(1,:));

%------------------------------------------
% Diffusion Process and new coordinate systemn

lamda_t=Lamda.^options.t;                    %Diffusion process by step t
% Z=lamda_t(2:dim_new+1)*(V(:,2:dim_new+1))';
Z=V(:,2:dim_new+1)*lamda_t(2:dim_new+1,2:dim_new+1);

%% Recording & Output
%-----------------------------------------
% Basic Parameter Information (Input Information)

% model.metric=options.metric;
% model.kernel=options.kernel;
% model.kpara=options.kpara;
% model.dim_new=dim_new;
% model.t=t;
model.options=options;
model.X=X;

%-----------------------------------------
% Processed Information (Result Information)
model.K=K;
model.D=D;              % Degree Matrix of K
model.L=L;
model.Ln=Ln;              % Degree Matrix of K
model.Lt=L^options.t;
model.dist=Distance;
model.eigenvalues=diag(Lamda);

model.V_dr=V(:,1:dim_new+1);                        %  Dimension Reduced Eigenvectors of L (The Laplacian Matrix)
model.Vn_dr=Vn(:,1:dim_new+1);                      %  Dimension Reduced Eigenvectors of Ln (The Normalized Laplacian Matrix)
model.Lamda_dr=Lamda(1:dim_new+1,1:dim_new+1);      %  Dimension Reduced Eigenvalue of L and Ln (They are the same)
model.C_dr=C(1:dim_new+1,1:dim_new+1);              %  Dimension Reduced Transformation Matrix between V and Vn

%-----------------------------------------
% Full Information. Conditional output. (For Further Research without recalculation. Mass Memoey needed)
if options.FullRec == 1
    model.V=V;              % Full Eigenvectors of L  (The Laplacian Matrix)
    model.Vn=Vn;            % Full Eigenvectors of Ln (The Normalized Laplacian Matrix)
    model.Lamda=Lamda;      % Full Eigenvalue of L and Ln (They are the same)
    model.C=C;              % Full Transformation Matrix between V and Vn
end

return


function [X_star,Anal,model] = Kpca_ReconAnal(X,KpcaOptions,PreiOprions)
% KPCA Kernel Principal Component Reconstruction and full Analysis. 
%  
% Synopsis:
% Synopsis:
% [X_star]=Kpca_ReconAnal(X,KpcaOptions,PreiOprions)
% [X_star]=Kpca_ReconAnal(X)
%
% Description:
% KPCA is run on X to achieve reduction dataset Z. Z is then used to 
% reconstruct X Kpca_PreImage.m function. Also this fuctions return all
% uesfull information to evaluate the performance of how where
% reconstruction is.
%
% Input:
%  X             [sample x dimension] indicates Training data.
%  KpcaOptions   [Structure]          Kpca parameters
%       .ker     [string] Kernel identifier (see 'help kernel'); 
%                (default 'linear').
%       .arg     [1 x narg] kernel argument; (default 1).
%       .new_dim [1x1] Output dimension (number of used principal 
%                components); (default num_data).
%  PreiOprions   [Structure]          Kpca PreImage solver parameter
%       .type      % Method of reconstruction form distance to position/
%                  % Choose from 'LSE'(Least square estimate) or
%                  % 'Dw'(Distance weight).                Default: LSE
%       .para      % Paratmter for Distance weight method. Default: 1
%       .neighbor  % Number of distances used. Choosing starts from
%                  % the shortest length.              Default: All
% Output:
%  Z [num x new_dim]            matrix of projections of X
%  model                        [structure] see below
%     .DR_method='Kpca';
%     .options = options;       [Structure]. Parameter for the Kpca of used model
%     .options.ker [string]     Kernel identifier (see 'help kernel').
%     .options.arg [1 x nargs]  Kernel argument(s).
%     .options.new_dim [1 x 1]  Output dimension (number of used principal
%     .X = X;                   [num_data x dim] original data 
%     .K=K;                     [num_data x num_data] Kernel matrix of X
%     .eigenvalues=Lambda;      all eigenvalue in vector form
%     .cputime = cputime - start_time;    % Time used for the process
%     .b [new_dim x 1]          Bias.
%     .Alpha [num x new_dim]    Multipliers.
%     
%      for FullRec
%     .eigenval=diag(Lambda);                           %Full U Eigenvectors     
%     .kercnt = num*(num+1)/2;
%     .MsErr = triu(ones(num,num),1)*model.eigval/num;  %MSE with respect to used basis vectors;
%     .mse = model.MsErr(options.new_dim);              %Mean square representation error of maped data   
%  Anal                         [structure] see below
%     .distance_star            [sample x sample] distance_star(i,j) means distance between ith Z_star and jth Z_star in original data speace. 
%  
% See also 
%  Kpca, Kernel, Kpca_PreImage.
%   
% About: 
%  Modification
%  WeiX, Dec 3rd 2014,First Edition
%
%
%% Initialization and Parameters
if nargin < 3, PreiOprions = []; end
if nargin < 2, KpcaOptions = []; end

if ~isfield(KpcaOptions,'ker'), KpcaOptions.ker='gaussian'; end             % Default output information;
if ~isfield(KpcaOptions,'arg'), KpcaOptions.arg=10000000; end               % Default output information;
if ~isfield(KpcaOptions,'new_dim'), KpcaOptions.new_dim=2; end              % Default output information;
if ~isfield(KpcaOptions,'FullRec'), KpcaOptions.FullRec = 0; end            % Default output information;

if ~isfield(PreiOprions,'type'), PreiOprions.type ='Dw'; end                % Default metric to measure the distance;
if ~isfield(PreiOprions,'para'), PreiOprions.para = 5; end                  % Default kernel function
% if ~isfield(PreiOprions,'neighbor'), PreiOprions.neighbor = num_data; end % Default using all

%% Main
[num,dim]=size(X);
X_star=zeros(num,dim);  %Assign memory

[Z_star,model] = Kpca(X,KpcaOptions);

%% Main      
for index=1:num %selecting every point in original space "X" to reconstruct  
    Z_starN=Z_star(index,:);  
    Do=zeros(num,1);      %Assign memory
    for i =1:num          % Distance_original space
        Do(i,1)  = Distance_OriginalSpace(Z_starN,i,model);    %Finding the distances original       
    end 
    X_star(index,:)= Dist2pos(model.X,Do,PreiOprions);
    Anal.Distance_star(index,:)=Do';
end       

Anal.Distance_orig = pdist(X,'euclidean');  % euclidean distance
Anal.Distance_orig = squareform(Anal.Distance_orig);
Anal.Distance_AbsDiff = abs(Anal.Distance_orig-Anal.Distance_star);
Anal.Distance_RateAbsDiff=Anal.Distance_AbsDiff./Anal.Distance_orig; %Avoide Inf and ignore the selfdistance defferences.
% Anal.Distance_RateAbsDiff=(Anal.Distance_AbsDiff-diag(diag(Anal.Distance_AbsDiff)))./(Anal.Distance_orig+eye(num)); %Avoide Inf and ignore the selfdistance defferences.


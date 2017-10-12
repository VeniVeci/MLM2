function [X_star]=Kpca_Recon(X,KpcaOptions,PreiOprions)
% Kpca Reconstruction via preimage solver
%
% Synopsis:
% [X_star]=Kpca_Recon(X,KpcaOptions,PreiOprions)
% [X_star]=Kpca_Recon(X)
%
% Description:
% KPCA is run on X to achieve reduction dataset Z. Z is then used to 
% reconstruct X Kpca_PreImage.m function.
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
% 
% Example:
%
% See also 
% Kpca, Kpca_PreImage
% 
% About: 
% 
% Modifications:
% WeiX, Dec-3rd-2014, first edition 

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

%% Reconstruction
[Z,model] = Kpca(X,KpcaOptions);
X_star = Kpca_PreImage(Z,model,PreiOprions);

        
return

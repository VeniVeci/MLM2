function [Z,model] = LppX(X,options)
% LppX: Locality Preserving Projections by XING
%
% Synopsis:
% [Z,model] = LppX(X)
% [Z,model] = LppX(X,option)
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
% Pcakage Require:
% Example:
%S
% About: 
%  Modification
%  WeiX, 4-11-2015, First Edition
%
%
%% Initialization and Parameters
start_time = cputime;
X=X';

[dim,num]=size(X);
if nargin < 2, options = []; end
if ~isfield(options,'neighborType'), options.neighborType = 'k'; end
if ~isfield(options,'neighborPara'), options.neighborPara = num/10; end % 10% of point as neighbor point;
if ~isfield(options,'metric'), options.metric ='euclidean'; end          % Default metric to measure the distance;
if ~isfield(options,'t'),options.tAuto = 1; end                             % Default parameter for heat kernel=3;
% if ~isfield(options,'kernel'), options.kernel = 'gaussian'; end       % Default kernel function
% if ~isfield(options,'kpara'),  options.kpara = 1000; end                 % Default kernel parameter
if ~isfield(options,'dim_new'),options.dim_new = 2; end                  % Default new dimension=3;
% if ~isfield(options,'t'), options.t = 1; end                             % Default Diffusion times;
if ~isfield(options,'FullRec'), options.FullRec = 0; end                 % Default output information;

%% Main
% ---------------------Construct the Graph---------------------------------
 DistM=pdist2(X',X',options.metric);     %distance matrix % Matlab statistic Tool box

% Auto choose distance regulator t

if options.tAuto==1
    options.t=sum(DistM(:).^2)/(num^2);
end
 
switch options.neighborType
    case 'k'
         num_neighbor=options.neighborPara;
         [D_E_sorted, D_E_sortIndex] = sort(DistM); 
         W_d=DistM;                     % Weight_distance(based)
         for i=1:num
             temp_coli=D_E_sorted(:,i);         % coli = column i th
             temp_coli(2+num_neighbor:end)=0;
             W_d(D_E_sortIndex(:,i),i)=temp_coli;
         end
         
         W_d=W_d';
         W_d = min(W_d,W_d');       %% Important. Ensure the matrix is symmetric
         
         clearvars temp_coli        % clean temp variable
         
    case 'epsilon'
        radius=options.neighborPara;
        W_d =  DistM.*(DistM<=radius); 
%       D = min(D,INF); 
    otherwise
        error('Error: Undefined type of neighbor.');    
 
end

G=(W_d~=0);
% W=X*X'.*G;
W=exp(-DistM.^2/options.t);     % Weight matrix
W=W.*G;
d=sum(W,2);                    
D=diag(d);                      % Degree Matrix  
L=D-W;                          % Laplacian Matrix

S_D = X * D * X';
S_L = X * L * X';
% S_D = (S_D + S_D') / 2;         % Ensure symetric
% S_L = (S_L + S_L') / 2;

% eigsoptions.disp = 0;
% eigsoptions.issym = 1;
% eigsoptions.isreal = 1;
% [eigvector, eigvalue] = eigs(LP, DP, no_dims, 'SA', eigsoptions);
[eigvector, eigvalue] = eig(S_L, S_D,'qz');


% [eigvalue, ind] = sort(diag(eigvalue), 'descend');
[eigvalue, ind] = sort(diag(eigvalue), 'ascend');

V=eigvector(:,ind(1:options.dim_new));
Z = V'*X;
Z=Z';


model.DR_method='LPP';
model.options = options;
% model.X = X;
model.V=V;      % Mapping basis

model.cputime = cputime - start_time;

% Full Information. Conditional output. (For Further Research without recalculation. Mass Memoey required)
if options.FullRec == 1
    model.eigenval=eigvalue;              
    model.eigenvec=eigvector;   
    model.cputime = cputime - start_time;     % Update process time
    model.G=G;
    model.W=W;
    
    
end

return


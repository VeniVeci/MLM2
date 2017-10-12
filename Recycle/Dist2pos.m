function X_star= Dist2pos(X,dist,options)
% Distances to positions coordinate calculation using advanced Different
% Schemes 
%
% Synopsis:
% X_star= Dist2pos(X,dist)
%
% Description:
% It computes the coordinate of a new point. The new ponit is added to 
% the original data space providing distances between the point and 
% each original data point.
% 
% 
% Input:
%       X [num_data X dim] is the original data 
%       dist [num_data x 1] is the distance between X_star and each X point  
%       options.
%               type  % Method of reconstruction form distance to position/
%                     % Choose from 'LSE'(Least square estimate) or
%                     % 'Dw'(Distance weight).                Default: LSE
%               para  % Paratmter for Distance weight method. Default: 1
%               neighbor  % Number of distances used. Choosing starts from
%                         % the shortest length.              Default: All
%
% Output:
%       X_star [1 x dim] Coordinate of X_star
%
% Example:
% See also:
% About: 
% Pcakage Require:
% See also:
% 
%
%  Modification
%  WeiX, Sep 13th 2014, First Edition
%
%% Initialization and Parameters
[num_data,dim]=size(X);
if nargin <= 2, options = []; end
if ~isfield(options,'type'), options.type ='LSE'; end          % Default metric to measure the distance;
if ~isfield(options,'para'), options.para = 1; end           % Default kernel function
if ~isfield(options,'neighbor'), options.neighbor = num_data; end           % Default kernel function

if options.neighbor ~= num_data
    [dist,index]=sort(dist);
    dist = dist(1:options.neighbor,:);
    X = X(index(1:options.neighbor),:);
end

[num_data,dim]=size(X);

%% Main 
switch options.type
    case 'LSE'                      % --LSE solution 
        dist=dist+(dist==0)*1e-50;  % Make sure the divided number is not zero. 1e-50 should be adjust according to machine.
        dist_square=dist.^2;    % Square Distance 
        means=mean(X);
        Xc=(eye(num_data,num_data)-ones(num_data,num_data)/num_data)*X; % X with centering
        
        do=sum(Xc.^2,2);
        d=do-dist_square;
        D=(d-(ones(num_data,num_data)/num_data*d))/2;
        %------------------------------------------
        % Choosing different Linear equation solver 
        %X_star = lsqlin(Xc,D);
        %X_star = lscov(Xc,D);
        X_star = Xc\D;              % A most common choice in MATLAB
        %------------------------------------------
        X_star=(X_star'+means);

    case 'Dw'                   % --Distance weight solution
        dist=dist+(dist==0)*1e-50;  % Make sure the divided number is not zero. 1e-50 should be adjust according to machine.
        k=options.para;         % Polynomial term of distance weight method
        %----------------------------------------------------------------------
        % Choosing different measure of similarity. Here is simple 1./(dist.^k)
        simi=1./(dist.^k);   
        simi_sum=sum(simi);
        weight=simi/simi_sum; % Weight of each point
        %----------------------------------------------------------------------
        weight=weight';
        X_star=weight*X;        % Weights times coordiantes
        
    otherwise
    error('Error: undefined distance to postion method.');

end

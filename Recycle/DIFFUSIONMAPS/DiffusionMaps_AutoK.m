function [ kpara ] = DiffusionMaps_AutoK(X,options)
% DiffusionMap Auto Kpara searching.
%
% Synopsis:
% [ kpara ] = DiffusionMap_AutoK(X,options)
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
%       kpara          % suitable kapra for  DiffusionMap.m              
%
% Pcakage Require:    DiffusionMap MK-II
% Example:
%
% See also 
% DiffusionMap_PreImage.m; DiffusionMap.m
% 
% About: 
%  Modification
%  WeiX, Nov 17th 2014, First edition
%  WeiX, Dec 4th  2014, Monir Update
%
%
%% Initialization and Parameters
% [num,dim]=size(X);
%if nargin < 2, options = []; else options=c2s(options); end
if ~isfield(options,'metric'), options.metric ='euclidean'; end          % Default metric to measure the distance;
if ~isfield(options,'kernel'), options.kernel = 'gaussian'; end    % Default kernel function
% if ~isfield(options,'kpara'),  options.kpara = 1000; end                 % Default kernel parameter
% if ~isfield(options,'dim_new'),options.dim_new = 3; end                  % Default new dimension=3;
% if ~isfield(options,'t'), options.t = 1; end                             % Default Diffusion times;
% if ~isfield(options,'FullRec'), options.FullRec = 0; end                 % Default output information;

%% Main
% Calculating distances
Distance =pdist2(X,X,options.metric);
if options.metric == 'euclidean'
   kpara=mean(Distance(:));
   kpara=kpara^2;
else
   err('Measurement Method not developed yet!')
end


end


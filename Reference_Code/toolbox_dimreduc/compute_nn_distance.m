function [D,nn_list] = compute_nn_distance(X,nbr_nn, options)

% compute_nn_distance - compute the distance to the nearest neighbors
%
%   [D,nn_list] = compute_nn_distance(X,nbr_nn, options);
%
% X is a (d,n) set of n points in R^d
% nbr_nn is the number of nearest neighbors to retrieve
%
% D is a (n,nbr_nn) matrix of distance, D(i,j) is distance between point
%   X(:,i) and point X(nn_list(i,j),:)
% nn_list(i,:) is the set of nearst neighbors
%
%   options.use_nntools = 0 force the use of the slow matlab code
%       for nearest neighbors computations.
%
%   options.exlude_self = 1 to avoid taking a point it self neighbor.
%
%   Copyright (c) 2006 Gabriel Peyré

options.null = 0;

if isfield(options, 'use_nntools')
    use_nntools = options.use_nntools;
else
    use_nntools = 1;
end
if isfield(options, 'exlude_self')
    exlude_self = options.exlude_self;
else
    exlude_self = 0;
end

if nargin<2
    nbr_nn = size(X,2);
end

if exlude_self==1
    nbr_nn = nbr_nn + 1;
end

nbr_nn = min(nbr_nn,size(X,2));

if exist('nn_prepare')>0 && use_nntools
    % use fast mex code
    atria = nn_prepare(X');
    [nn_list,D] = nn_search(X', atria, X', nbr_nn, 0);
else
    % use slow matlab code
    D1 = sqrt( compute_distance_matrix(X) );
    % find closest points
    [D,nn_list] = sort(D1);
    D = D(1:nbr_nn,:)';
    nn_list = nn_list(1:nbr_nn,:)';
end

if exlude_self==1
	% remove self reference
	nn_list = nn_list(:,2:end);
	D = D(:,2:end);
end
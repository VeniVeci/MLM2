function mappedX = mds(X, no_dims)
%MDS Run MDS on the data to get a low-dimensional visualization
% 
%   mappedX = mds(X, no_dims)
%
% Run multidimensional scaling on the dataset X to get a two-dimensional 
% visualization. The low-dimensional representation is returned in mappedX.
% It has dimensionality no_dims (default = 2).
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.6b.
% The toolbox can be obtained from http://www.cs.unimaas.nl/l.vandermaaten
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten
% Maastricht University, 2008


    if ~exist('no_dims', 'var')
        no_dims = 2;
    end

    % Initialize some variables
	mappedX = compute_mapping(X, 'PCA', no_dims);
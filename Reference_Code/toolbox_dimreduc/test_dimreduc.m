% simple tests for dimension reduction on 3D data
%
%   Copyright (c) 2005 Gabriel Peyré


name = 'square';
name = '3d_cluster';
name = 'swissroll';
name = 'puncted_sphere';
name = 'scurve';


n = 1000;
options.sampling = 'uniform';
options.sampling = 'rand';
options.noise_level = 0;
[X,col] = load_points_set( name, n, options );

use_lle = 1;
use_isomap = 1;
use_hlle = 1;
use_leigs = 1;
use_ltsa = 1;
use_pca = 1;

% number of neighbors
nn_nbr = 12;

xy = {}; lgd = {};
xy{1} = X;
lgd{1} = '';

if use_isomap
    disp('----- ISOMAP -----');
    clear options;
    options.nn_nbr = nn_nbr;
    xy{end+1} = isomap(X,2, options);
    lgd{end+1} = 'Isomap';
end
if use_lle
    disp('----- LLE -----');
    nn_nbr = nn_nbr;
    xy{end+1} = lle(X,nn_nbr,2);
    lgd{end+1} = 'LLE';
end
if use_hlle
    disp('----- HLLE -----');
    xy{end+1} = hlle(X,nn_nbr,2);
    lgd{end+1} = 'HLLE';
end
if use_leigs
    disp('----- Laplacian Eigenmaps -----');
    xy{end+1} = leigs(X, 2, nn_nbr);
    lgd{end+1} = 'Laplacian Eigenmaps';
end
if use_ltsa
    disp('----- LTSA -----');
    xy{end+1} = ltsa(X,2,nn_nbr);
    lgd{end+1} = 'LTSA';
end

if use_pca
    disp('----- PCA -----');
    [Y,xy{end+1}] = pca(X,2)
    lgd{end+1} = 'PCA';
end

%%% display %%%
p = length(xy);
nb = ceil(sqrt(p));
clf;
for i=1:p
    subplot(nb,nb,i);
    plot_scattered(xy{i},col);
    axis off;
    title(lgd{i});
end
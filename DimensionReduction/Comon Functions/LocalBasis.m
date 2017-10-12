function [ V ] = LocalBasis( X, Z )
% Local basis finder try to find a local basis describe connection between
% X and Z. The basis are NOT orthonormal 
%
% Report: it is not working well. Unstable. Need more insvigate.
% 
% Input:
% X                         % [dimension_x x sample] original Training data.
% Z                         % [dimension_z x sample coefficients/representations 
%                           for PCA. dim_new <= dimension !!
% Output:
% V                         % [dimension_x x number of basis]  Basis
%
%
% Modifications:
% WeiX, 27-10-2015, first edition. Conclusion:Fail.


% V*Z'=X

V = (Z'\X')';




end


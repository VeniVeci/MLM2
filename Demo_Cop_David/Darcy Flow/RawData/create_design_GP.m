
% 25 January 2017

% This routine uses precomputed eigenmodes (of the correlation function) for creating a 
% training set for the GP emulator 

% Load precomputed data
load ./RawData/Eigenmodes/lambda_0_3/S_elliptic.mat;      % precomputed eigenvalues
load ./RawData/Eigenmodes/lambda_0_3/A_elliptic.mat;      % precomputed eigenvectors

% ----------------------    User specifications    ----------------------%
GPvariance = 1.0;              % variance in the correlation formula
GPmean = 0.0;                  % mean of Z, log K = Z
m=8;
ndesign = 2^m; %=256;          % Number of design points (use power os 2 for Sobol sequences)
ntrunc =  length(S(:,1));      % dimension of the basis (number of eigenvectors)

kl_trunc = 1111;               % 97.91% variance preservation             
                               % number of kl coefficients to be retained, there's a current 
                               % limitation in Matlab. You can use upto 1111 dimensions. 
                               % The NAG toolbox allows to generate more but I do not have access 
                               % to it anymore. For now, let's try with the
                               % maximum allowed then. 
% ---------------------------------------------------------------------- %

Nz=1;Dx=1;Dy=1;Dz=1;           % variables related to the computational domain

% Sample points from a Sobol sequence
sobolvar = sobolset(kl_trunc); 
design = net(sobolvar, ndesign+1)'; 

% Remove the point at zero
design = design(:, 2:(ndesign+1));

% Push sampled points forward a N(O,I) distribution
for i=1:kl_trunc
    design(i,:) = norminv(design(i,:), 0, 1.0); % NO BROADER Design, for now use N(0,I) design points
end

% compute log K = Z
Zdesign = sqrt(GPvariance) * A(:,1:kl_trunc) * sqrt(S(1:kl_trunc,1:kl_trunc)) * design + GPmean;    
% Compute K
Kdesign = exp(Zdesign);       

    %Modified By WeiX to test
    for i=1:ndesign
        iZdesign=Zdesign(:,i);
        iZdesignMatrix=reshape(iZdesign,50,50);
        iKdesignMatrix = exp(iZdesignMatrix);  
        Kdesign2(:,i)=iKdesignMatrix(:);
    end




save('./RawData/data_lambda_0_3_sigma_1_0/design','design');     % design

save('./RawData/data_lambda_0_3_sigma_1_0/K_design','Kdesign');  % Matrix of Nl design 
                                                                    % permeabilities
                                                                    % (columns)
                                                                    
% ***********************************  END CODE ******************************************* %                                                                   
     





% Last modified: 25 Jan 2017 

% This routine generate Nl permeability fields from the precomputed eigen
% decomposition of the correlation matrix C, i.e., matrices A and S such that C=A*S*A^(-1)  


% User specifications
Nl = 1000;                  % Number of TRUE permeability fields to generate
ntrunc = 2499;              % Number of eigenmodes used in precomputed A and S
kl_trunc = ntrunc;          % Number of KL coefficients retained in the KL expansion

GPvariance =  1.0;          % variance of the random fields (gpvar)              
GPmean = 0.0;
%GPmean = -GPvariance./2;   % This formula will lead to a log normal
                            % distribution of mean 1

% load precomputed matrices A and S
load('./RawData/Eigenmodes/lambda_0_3/A_elliptic.mat'); % eigenvectors
load('./RawData/Eigenmodes/lambda_0_3/S_elliptic.mat'); % eigenvalues

u = randn(kl_trunc,Nl);   % KL coefficients    \xi ~ N(0,1)                                                                        
Z = sqrt(GPvariance) * A(:,1:kl_trunc) * sqrt(S(1:kl_trunc,1:kl_trunc)) * u + GPmean;        
K = exp(Z); 

save('./RawData/data_lambda_0_3_sigma_1_0/kl_coeff','u');        % Matrix of Nl permeabilities
save('./RawData/data_lambda_0_3_sigma_1_0/K_true','K');  % Matrix of Nl permeabilities

% -----------------------   Plot permeability  ---------------------------
perme_num = 5; 

Nx=50;Ny=50;    % nodes each direction
Dx=1;Dy=1;      % [0,1]x[0,1]

% Check the min and max for other perme_fields for   caxis([ ]);
% min(min(perme_field)) 
% max(max(perme_field)) 

% Plot Permeability field
Kmain=reshape(K(:,perme_num),50,50);

%figure(perme_num);
figure(1);
pcolor([0:Dx/(Nx-1):Dx],[0:Dy/(Ny-1):Dy],Kmain);
shading interp;
%shading flat;
%caxis([0.0 10.0]); % min(min(perme_field))  and max(max(perme_field))
colormap jet;
colorbar('FontSize',14);
xlabel('$x$','FontSize',18,'interpreter','latex');
ylabel('$z$','FontSize',18,'interpreter','latex');


% *************************** END CODE ***********************************

























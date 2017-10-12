
% 25 January 2017

% This routine performs the eigen decomposition of the correlation function  

% ------------------ Specification of Mesh and Domain --------------------
Nx=50;    % Number of Nodes axis-x
Ny=50;    % Number of Nodes axis-y
Nz=1;     % 2-D plane
Dx=1;   % length of domain
Dy=1;   % height of domain
Dz=1;     % 2-D plane
% ------------------ Gaussian Process Specifications ---------------------
% GPvariance =  variance of the random fields (gpvar) is specified in the KL expansion
              
% GPmean = -GPvariance./2;     % Same above, this formula will lead to a log normal
                               % distribution of mean 1. For now we are setting GPmean=0

lenscale = 0.3;


ntrunc =  2499;                % number of KL basis vectors to use
                               % this choice guarantees the 100% of the total variance
               
% -------------------- CORRELATION FUNCTION -----------------------------

%  Nodes where to evaluate Permeability at 
[X Y]=meshgrid(0.02:0.02:Dx, 0.02:0.02:Dy);

X=reshape(X,Nx*Ny,1);
Y=reshape(Y,Nx*Ny,1);

x=[Y X];

n = size(x,1);
     
% Correlation function 
k = @(x,y) exp(-sqrt((x-y)*(x-y).')/lenscale);       

C = zeros(n,n);
for i = 1:n
    for j = 1:n
        C(i,j) = k(x(i,:), x(j,:));
    end
end

[A,S] = eigs(C,ntrunc);   % eigen decomposition of the correlation matrix

save('./ZZZ_Third_paper/Eigenmodes/lambda_0_3/S_elliptic','S'); % eigenvalues
save('./ZZZ_Third_paper/Eigenmodes/lambda_0_3/A_elliptic','A'); % eigenvectors

% Truncation error
trunc_error = sum(diag(S))./(Nx*Ny) % truc_error must be close to 1

% ***********************************  END CODE ******************************************* % 































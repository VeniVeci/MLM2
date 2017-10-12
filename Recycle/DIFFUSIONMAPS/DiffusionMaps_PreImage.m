function [X_star] = DiffusionMaps_PreImage(Z_star,model,options)
% function of Diffusion dimension reduction preimage solution.      MARK-II
%
% Synopsis:
% X_star = DiffusionMap_PreImage(Z_star,model). Default 'LSE' distance to coordinate method.
% X_star = DiffusionMap_PreImage(Z_star,model,options). 
% [X_star,Dist_star] = DiffusionMap_PreImage(Z_star,model,options) Output with recover distance informations.
%
% Description:
% The function find a new point's position in original dataspace using the
% offered position information in embedded space.
% 
% steps for the program:
% 
% Input:
% Z_star [Samples X Dimensions]         The new point in embedded space
% model [structure]
%      .options=options;        % Diffusion map parameters
%      .X=X;                    % Original dataset [Samples X Dimensions]
%      .K=K;                    % Kernel Matrix    [Samples X Samples]
%      .D=D;                    % Degree Matrix of K
%      .L=L;                    % The Laplacian Matrix
%      .Ln=Ln;                  % The Normalized Laplacian Matrix
%      .Lt=L^options.t;         % Laplacian Matrix after t time.
% 
%      .V_dr=V(:,1:dim_new+1);                        %  Dimension Reduced Eigenvectors of L (The Laplacian Matrix)
%      .Vn_dr=Vn(:,1:dim_new+1);                      %  Dimension Reduced Eigenvectors of Ln (The Normalized Laplacian Matrix)
%      .Lamda_dr=Lamda(1:dim_new+1,1:dim_new+1);      %  Dimension Reduced Eigenvalue of L and Ln (They are the same)
%      .C_dr=C(1:dim_new+1,1:dim_new+1);              %  Dimension Reduced Transformation Matrix between V and Vn
% 
% PreIoptions [structure]           % Options for the Pre-image solution      
%        .type                 % Type of distance to coordinate method.   'Dw'(Distance weight)/'LSE'(Least square estimate
%        .para                 % Parameter of distance to coordinate recover method. For 'Dw' method Only.
% 
%
% Output:
% X_star [samples x dimension_new]  % Coordinate of X_star in original data
%                                     space, corresponding to Z_star.
%
% Pcakage Require: DiffusionMap MK-II
% Example:
%
% See also 
% DiffusionMaps.m;
% 
% About: 
% 
% Modifications:
%  WeiX, Sep 13th 2014, First Edition
%  WeiX, Dec 4th  2014, Minor Update
%% ----------------------Initialization------------------------------------
if nargin <= 2, options = []; end


[Num_Zstar,Dim_Zstar]=size(Z_star);
[Num_X,Dim_X]=size(model.X);

X_star=zeros(Num_Zstar,Dim_X);          % Assign memory
Dist_star=zeros(Num_Zstar,Num_X);       % Assign memory

Lamda=model.Lamda_dr;
C=model.C_dr;
D=model.D;
V=model.V_dr;    
v11=model.V_dr(1,1);                     % First eigenvector value. They are all the same.

switch model.options.kernel
    case 'gaussian'
        k_starstar = exp(-(0).^2/model.options.kpara);    % Actually the "0" could change due to differnet metric method.       
    otherwise
        error('Error: Undefined type of kernel function.');    
end

%% -----------------------Main--------------------------------
for i = 1:Num_Zstar

%     Z_star_complete=[v11*Lamda(1,1),Z_star(i,:)];  
% %     v_star=Z_star(i,:)*Lamda(1:end,1:end)^(-model.options.t);
%     v_star=Z_star_complete*Lamda(1:end,1:end)^(-model.options.t);
    
    v_starM1=Z_star(i,:)*Lamda(2:end,2:end)^(-model.options.t);  %V_starM1= Vectorc_star minus the ist component
    v_star=[v11,v_starM1];
    
    d_star=sqrt(k_starstar)/sqrt(v_star*C^2*Lamda*v_star'); 
    % k_starstar=d_star*v_star*C^2*Lamda*v_star'*d_star;
    k_star=D*V*C^2*Lamda*v_star'*d_star;
    
%     K_orig=D*V*C^2*Lamda*V'*D;

% -------------!!!!!!!!! Correction !!!!!!!!!!!!-------------------------- 
k_star=(k_star>=0).*k_star;             % Ensure k_star is positive.
k_star=k_star+(k_star==0).*1e-50;       % Ensure 0 in k_star are limt to 0 Rather than real '0'..
% ------------------------------------------------------------------------

switch model.options.kernel
    case 'gaussian' 
        dist_star= sqrt(-log(k_star)*model.options.kpara);
    otherwise
        error('Error: Undefined type of kernel function.');    
end

%     if nargout > 1
%         dist_star = dist_star;        
%     else

%    Dist_star(i,:)=dist_star';     Record the distance for further use. Usually for further exam.
    
    X_star(i,:)= Dist2pos(model.X,dist_star,options);
    
end
    


return % End of Function



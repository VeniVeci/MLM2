function Distance_orig = Distance_OriginalSpace(Z_star,i,model)
% Distance computation between two ponit in origial data space
%  
% Synopsis:
% Distance_orig = Distance_OriginalSpace(Z_star,i,model)
%
% Description:
%  It computes distance between vectors mapped into the feature 
%  space induced by the kernel function (model.options.ker,
%  model.options.arg).The distance is computed between projection
%  mapped from feature space and the ith point in the feature space 
%  given by model:
% 
% steps for the program:
%   d_ij^2=log((1/2)*(K_ii+K_jj-D_ij^2)) *(-2*C^2)
%   In gaussian kernel the K_ii=K_jj=0
%   d_ij is the distance in original data space between point i and j
%   D_ij is the distance in feature space between point i and j
%   C is the parameter of gaussian kernel
% 
% Input:
%  Z_star [Dim x 1] 
%  i [1,1] ith observation of X
%  model                        [structure] see below
%     .DR_method='Kpca';
%     .options = options;       Structure. Parameter for the Kpca of used model
%     .options.ker [string]     Kernel identifier (see 'help kernel').
%     .options.arg [1 x nargs]  Kernel argument(s).
%     .options.new_dim [1 x 1]  Output dimension (number of used principal
%     .X = X;                   [num_data x dim] original data 
%     .K=K;                     [num_data x num_data] Kernel matrix of X
%     .eigenvalues=Lambda;      all eigenvalue in vector form
%     .cputime = cputime - start_time;    % Time used for the process
%     .b [new_dim x 1]          Bias.
%     .Alpha [num x new_dim]    Multipliers.
%     
%      for FullRec
%     .eigenval=diag(Lambda);                           %Full U Eigenvectors     
%     .kercnt = num*(num+1)/2;
%     .MsErr = triu(ones(num,num),1)*model.eigval/num;  %MSE with respect to used basis vectors;
%     .mse = model.MsErr(options.new_dim);              %Mean square representation error of maped data.
% 
% Output:
% Distance [1 x 1] distance calculated
% 
% Example:
%
% See also 
% Kernel, Kpca. 
% 
% About: 
% 
% Modifications:
% 19-jul-2013, WeiX, first edition 
% WeiX, Dec 1st  2014, MK-II
% WeiX, 13-11-2015, Add ability to handle Dimension-dismatch Z_star

[num_data,dim_X] = size(model.X);
[num_Z,dim_Z] = size(Z_star);    
% if dim_Z~=model.options.new_dim
%     warning('Dimension dismatch. Program automatically truncate to match');
% end


if strcmp(model.options.ker,'linear')

elseif strcmp(model.options.ker,'poly') 

    [num_data,dim_X] = size(model.X);
    Beta_star=model.Alpha(:,1:dim_Z)*Z_star'+ones(num_data,1)/num_data;
    Kyi=model.K(i,:)';

    K_istar=Beta_star'*Kyi;
    K_ii=model.K(i,i);
    K_starstar=Beta_star'*model.K*Beta_star;

    d=model.options.arg(1);
    c=model.options.arg(2);

    XiXstar = K_istar^(1/d)-c;
    XiXi = K_ii^(1/d)-c;
    XstarXstar = K_starstar^(1/d)-c;

    Distance_orig_square=XiXi+XstarXstar-2*XiXstar;
    
    Distance_orig=sqrt(Distance_orig_square);
    
elseif strcmp(model.options.ker,'gaussian')   
    
    [num_data,dim_X] = size(model.X);
%     Beta=model.Alpha*Z_star'+ones(num_data,1)/num_data;
    Beta=model.Alpha(:,1:dim_Z)*Z_star'+ones(num_data,1)/num_data;
    Kyi=model.K(i,:)';
    Distance_feature_square=model.K(i,i)+Beta'*model.K*Beta-2*Beta'*Kyi;
    
    C=model.options.arg(1);
    Distance_orig_square=log((1/2)*(2-Distance_feature_square))*(-2*C^2);
    if Distance_orig_square<0 
        warning('Distance square in original space is minus number.');
    end
    if ~isreal(Distance_orig_square)
%         Distance_orig_square=real(Distance_orig_square);
        warning('Distance in original space is non-real number.');
    end
    
   
    Distance_orig=sqrt(Distance_orig_square);
%     Distance_orig=realsqrt(Distance_orig_square);
    
    
   
 elseif strcmp(model.options.ker,'sigmoid')  
     
    [num_data,dim_X] = size(model.X);
%     Beta_star=model.Alpha*Z_star'+ones(num_data,1)/num_data;
    Beta_star=model.Alpha(:,1:dim_Z)*Z_star'+ones(num_data,1)/num_data;
    Kyi=model.K(i,:)';

    K_istar=Beta_star'*Kyi;
    K_ii=model.K(i,i);
    K_starstar=Beta_star'*model.K*Beta_star;

    alpha=model.options.arg(1);
    c=model.options.arg(2);
   
    XiXstar = (atanh(K_istar)-c)/alpha;
    XiXi = (atanh(K_ii)-c)/alpha;
    XstarXstar = (atanh(K_starstar)-c)/alpha;

    Distance_orig_square=XiXi+XstarXstar-2*XiXstar;
    
    Distance_orig=sqrt(Distance_orig_square);
     
     
else
    disp('Error: undefined kernel.');
end

return
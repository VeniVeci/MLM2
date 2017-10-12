function X_star = Kpca_PreImage(Z_star,model,options)
% Kpca_Preimage solver.
%
% Synopsis:
%  X_star = Kpca_PreImage(Z_star,model)
%
% Description:
% Solve the Kpca preimage problem by a linear method introduced by JTY
% Kwok.
% It computes distance between vectors mapped into the feature 
% space induced by the kernel function (model.options.ker,
% model.options.arg).The distance is computed between projection
% mapped from feature space and the ith point in the feature space 
% given by model:
%
% Input:
%  Z_star [num x new_dim]       the point in reduced space without known original coordinate
%     model                     [structure] see below
%     .DR_method='Kpca';
%     .options = options;       Structure. Parameter for the Kpca of used model [usefuf information here]
%     .options.ker [string]     Kernel identifier (see 'help kernel').
%     .options.arg [1 x nargs]  Kernel argument(s).
%     .options.new_dim [1 x 1]  Output dimension (number of used principal
%     .X = X;                   [num_data x dim] original data [usefuf information here]
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
% options
%     .type      % Method of reconstruction form distance to position/
%                % Choose from 'LSE'(Least square estimate) or
%                % 'Dw'(Distance weight).                Default: LSE
%     .para      % Paratmter for Distance weight method. Default: 1
%     .neighbor  % Number of distances used. Choosing starts from
%                % the shortest length.              Default: All
%
% Output:
%  X_star [num x dim] corresponding coordinates in original space of Z_star
%
% Example:
% 
% See also 
% Kpca, kernel
% 
% Modifications:
% 19-jul-2013, WeiX, first edition 
% WeiX, Dec 1st  2014, MK-II
%
% Reference:
% The pre-image problem in kernel methods.(paper)
% Statistical Pattern Recognition Toolbox 1999-2003, Written by Vojtech Franc and Vaclav Hlavac

%% Initialization and Parameters
% timer

[Num_Zstar,Dim_Zstar]=size(Z_star);
[Num_X,Dim_X,]=size(model.X);
X_star=zeros(Num_Zstar,Dim_X);  %Assign memory

%% Main      
for index=1:Num_Zstar %selecting every point in original space "X" to reconstruct 
 
    Z_starN=Z_star(index,:);  
    Do=zeros(Num_X,1);      % Assign memory
    for i =1:Num_X          % Distance_original space
        Do(i,1)  = Distance_OriginalSpace(Z_starN,i,model);    %Finding the distances        
    end 
    
    % -------------!!!!!!!!! Correction !!!!!!!!!!!!-------------------------- 
    if ~isreal(Do)
        Do=real(Do);
        warning('Distance in original space is non-real number. Forced correction uisng the real part');
    end
    
%     Do=(Do>=0).*Do;             % Ensure k_star is positive.
%     Do=Do+(Do==0).*1e-50;       % Ensure 0 in k_star are limt to 0 Rather than real '0'..
    %----------------------------------------------------------
  
    switch options.type
        case 'LpcaI'
%             options.InMethod = 'LSE';          		 % Default Interpretation method;  
%             options.dim_new=2;
            X_star(index,:)= Pos2PosLpcaI(model.X,model.Z,Do',Z_starN,options);
        otherwise
            X_star(index,:)= Dist2pos(model.X,Do,options);
    end
    
end       

return



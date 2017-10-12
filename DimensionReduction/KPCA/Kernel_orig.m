function K = Kernel_orig(X,kernel,arg)
% KERNEL Evaluates kernel function.
%
% Synopsis:
%  K = kernel(X,kernel,arg)
%
% Description:
%  K = kernel( X, ker, arg ) returns kernel matrix K [n x n] 
%
%    K(i,j) = k(X(:,i),X(:,j))  for all i=1..n, j=1..n,
%     
%   kernel        Name           Definition
%   'linear'  ... linear kernel  k(a,b) = a'*b
%   'poly'    ... polynomial     k(a,b) = (a'*b+arg[2])^arg[1]
%   'rbf'     ... RBF (Gaussian) k(a,b) = exp(-0.5*||a-b||^2/arg[1]^2)
%   'sigmoid' ... Sigmoidal      k(a,b) = tanh(arg[1]*(a'*b)+arg[2])
%
% Input:
%  X [dim x n] Single matrix of input vectors.
%  options [struct] Decribes kernel and output dimension:
%   .ker [string] Kernel identifier (see 'help kernel'); 
%     (default 'linear').
%   .arg [1 x narg] kernel argument; (default 1).
%
% Output:
%  K [n1 x n1] or K [n1 x n2] Kernel matrix.
% 
% Example:
%  X = rand(2,50);
%  K = kernel( X, 'rbf', 1);
%  figure; pcolor( K );
%
% See also:
% Kpca
%
% Reference:
% Statistical Pattern Recognition Toolbox 1999-2003, Written by Vojtech Franc and Vaclav Hlavac

[dim,num_data] = size(X);

K=zeros(num_data,num_data);

if strcmp(kernel,'linear')
    K = X'*X;
%    for i=1:num_data
%        for j=1:num_data
%            K(i,j) = X(:,i)'*X(:,j);
%        end
%    end
    
elseif strcmp(kernel,'poly') 
    for i=1:num_data
        for j=1:num_data
            K(i,j) = (X(:,i)'*X(:,j)+arg(2))^arg(1);
        end
    end
    
%     K = (X'*X+ones(num_data)*arg(2)).^arg(1);
    
    
elseif strcmp(kernel,'rbf')   
    
    % An advance algorihm to calculate K 3 times faster compared with
    % regular one. !!! But it costs num_data times space
    %{
    for i=1:num_data
        tmp(:,i) = sum(bsxfun(@minus,X,X(:,i)).^2,1);
        K(:,i) = exp(-0.5*tmp(:,i)/arg(1)^2);
    end    
    %}
    
    %regular algorihm to calculate K
    for i=1:num_data
        for j=1:num_data
            temp=sum((X(:,i)-X(:,j)).^2);
            %K(i,j)=sum(X(:,i).^2+X(:,j).^2-2*X(:,i).*X(:,j));
            %K(i,j)=norm(X(:,i)-X(:,j)).^2;         
            K(i,j) = exp(-0.5*temp/arg(1)^2);
            
            %sum(X(:,i).^2+X(:,j).^2-2*X(:,i).*X(:,j));
            % sum((X(:,i)-X(:,j)).^2);
            %k(a,b) = exp(-0.5*||a-b||^2/arg[1]^2)
        end
    end
    
    % A better version ------------
%     D = pdist2(X',X');
%     K = exp((D.^2)./(-2*para(1)^2));
    % -------------------------------
    
 elseif strcmp(kernel,'sigmoid')  
    for i=1:num_data
        for j=1:num_data            
            K(i,j) = tanh(arg(1)*(X(:,i)'*X(:,j))+arg(2));        
        end
    end   
else
    disp('Error: undefined kernel.');
end

return

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 19-sep-2004, VF
% 5-may-2004, VF

% MEX-File function.                        
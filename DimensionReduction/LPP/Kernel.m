function K = Kernel(X,kernel,para)
% KERNEL Evaluates kernel function.
%
% Synopsis:
% K = Kernel(X,kernel,para)
% K = Kernel(X)
%
% Description:
% returns kernel matrix K [n x n]  where K(i,j) = k(X(:,i),X(:,j))  for all i=1..n, j=1..n,
%     
%   kernel        Name           Definition
%   'linear'  ... linear kernel  k(a,b) = a'*b
%   'poly'    ... polynomial     k(a,b) = (a'*b+para[2])^para[1]
%   'gaussian'... Gaussian       k(a,b) = exp(-0.5*||a-b||^2/para[1]^2)
%   'sigmoid' ... Sigmoidal      k(a,b) = tanh(para[1]*(a'*b)+para[2])
%
% Input:
%  X      [n x dim]          Single matrix of input vectors.
%  kernel [string]           Kernel identifier (see 'help kernel'); (default 'linear').
%  para   [1 x n*parameter]  kernel argument; (default 1).
%
% Output:
%  K [n x n] Kernel matrix.
% 
% Example:
%  X = rand(2,50);
%  K = kernel( X', 'rbf', 1);
%  figure; pcolor( K );
%
% See also:
% Kpca
%
%  Modification
%  WeiX, Lost date   , First Edition
%  WeiX, Dec 1st 2014, Improved algorithm.
%
% Reference:
% Statistical Pattern Recognition Toolbox 1999-2003, Written by Vojtech Franc and Vaclav Hlavac

%% Initialization and Parameters
[num,dim]=size(X);

if nargin <= 2
    para=[1,1];
end
if nargin <= 1
    kernel='linear'; 
end

K=ones(num,num);

switch kernel
    case 'linear'
        K = X*X';
        
    case 'poly'
%         % A slow version
%         for i=1:num
%             for j=1:num
%                 K(i,j) = (X(:,i)'*X(:,j)+para(2))^para(1);
%             end
%         end        

        % An improved version
        K = (X*X'+para(2)).^para(1);
        
    case 'gaussian'
        % 1)An advance algorihm to calculate K 3 times faster compared with
        % regular one. !!! But it costs num_data times space
        %{
        for i=1:num_data
            tmp(:,i) = sum(bsxfun(@minus,X,X(:,i)).^2,1);
            K(:,i) = exp(-0.5*tmp(:,i)/arg(1)^2);
        end    
        %}

        % 2) regular algorihm to calculate K
%         for i=1:num
%             for j=1:num
%                 temp=sum((X(:,i)-X(:,j)).^2);
%                 %K(i,j)=sum(X(:,i).^2+X(:,j).^2-2*X(:,i).*X(:,j));
%                 %K(i,j)=norm(X(:,i)-X(:,j)).^2;         
%                 K(i,j) = exp(-0.5*temp/para(1)^2);
% 
%                 %sum(X(:,i).^2+X(:,j).^2-2*X(:,i).*X(:,j));
%                 % sum((X(:,i)-X(:,j)).^2);
%                 %k(a,b) = exp(-0.5*||a-b||^2/arg[1]^2)
%             end
%         end
        
        % 3) a more efficient way in matlab
        D = pdist2(X,X);
        K = exp((D.^2)./(-2*para(1)^2));

        
        
    case 'sigmoid'
          % slow version
%         for i=1:num
%             for j=1:num            
%                 K(i,j) = tanh(para(1)*(X(:,i)'*X(:,j))+para(2));        
%             end
%         end   


        % An improved version
        K = tanh(X*X'*para(1)+para(2));
        
        
    otherwise
        disp('Error: undefined kernel.');

end

end
       

















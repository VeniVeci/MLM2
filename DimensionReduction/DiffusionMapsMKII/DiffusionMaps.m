function [Z,model] = DiffusionMaps(X,options)
% Description:
% Diffusion map Dimension Reduction Function. Version MK-II
%
% Synopsis:
% [Z,model] = Diffusionmap(X)
% [Z,model] = Diffusionmap(X,option)
%
% Input: 
%       X              % Original dataset [samples X dimension ] matrix
%       option         % Options of the diffusion map process
%             .kernel  % Type of kernel function
%             .kpara   % Parameter of the kernel function
%             .metric  % Method of measurement. Metric
%             .t       % optional time parameter in diffusion map 
%                          (default: model with multiscale geometry)    
%             .FullRec % Record all output information flag. 0 for false
%                          and 1 for ture.
%             .Ztype   % Output Z type. With (1) or without (0) 1st
%                           component (Theoritically constact)
%
% Output: 
%       Z              % new coordinates system [new dimension x samples] 
%       model          % preserving the Diffusion map process information and
%                      parameters
%
% About: 
%  Modification
%  WeiX, 11-11-2015, First Edition
%  WeiX, 4-1-2016,  Add Type2 to show Z

%
%
%% Initialization and Parameters
% [num,dim]=size(X);
if nargin < 2, options = []; end
if ~isfield(options,'metric'), options.metric ='euclidean'; end          % Default metric to measure the distance;
if ~isfield(options,'kernel'), options.kernel = 'gaussian'; end          % Default kernel function
if ~isfield(options,'kAuto'),  options.kAuto=0; end                      % Default no auto para
if ~isfield(options,'kpara'),  options.kAuto=1; end                      % Default kernel parameter auto searching
if ~isfield(options,'dim_new'),options.dim_new = 2; end                  % Default new dimension=3;
if ~isfield(options,'t'), options.t = 1; end                             % Default Diffusion times;
if ~isfield(options,'FullRec'), options.FullRec = 0; end                 % Default output information;
if ~isfield(options,'Ztype'), options.Ztype = 1; end                     % 0 without 1st component/ 1 with 1st component

dim_new=options.dim_new;
[num,dim]=size(X);

%% Main
% Calculating distances
Distance =pdist2(X,X,options.metric);

% Calculating the Kernel Matrix
switch options.kernel
    case 'gaussian'       
        if options.kAuto==1
            options.kpara=sum(Distance(:).^2)/(num^2);
        end
        K = exp(-Distance.^2/options.kpara);            
    otherwise
        error('Error: Undefined type of kernel function.');    
end

d=sum(K,2);                             % row sum
D=diag(d);                              % Degree Matrix. It is a diagonal matrix. 
P=D^(-1)*K;                             % The Markov matrix
Pprime=D^(-0.5)*K*D^(-0.5);             % The Normalized Markov matrix

% Eigen Decomposition of The Normalized  Markov matrix
[S,LamdaS] = eig(Pprime);  
[LamdaS,Index] = sort(diag(LamdaS),'descend');
S = S(:, Index);
LamdaS=diag(LamdaS);                     % LamdaS must = Lamda

R=D^(-0.5)*S;  % R2 must =R     % Right eigenvalue of P
L=D^(0.5)*S;                    % Left  eigenvalue of P


% %OLD VERSION------------------------------
%     % Check point! Not necessary for algorithm---------------------------------
%     % Eigen Decomposition of The Markov matrix
%     [R2,Lamda] = eig(P);  
%     [Lamda,Index] = sort(diag(Lamda),'descend'); % sort eigenvalues in descending order
%     R2 = R2(:, Index);                             % sort eigenvector in corresponding order to its eigenvalue
%     Lamda=diag(Lamda);
% 
% 
%     % Detection.---------------------------------------------------------------
%     % This part is made to ensure DiffusionMap is carried out safely
%     % under the limit of computer. Situation happens when K(i,j)=0.
%     if abs(mean(R2(:,1))-R2(1,1))>1e-3  % then each element in V are assumed to be the same 
%         error('parameter "options.kpara" in Diffusion Map is not suitable(Probably too small). Please Change.') 
%     end
% 
%     if ~isreal(R2)
%        warning('Eigenvector V of Laplacian Matrix L is not real number. Result might not be accurate. Parameter "options.kpara" is adviced to change. This could be normal as the lst few value of V tend to be complex number vector') 
%     end
%     % if sum(K(:,1))>length(K(:,1))/3
%     %     warning('Kernel Matrix varies very limit. Parameter "options.kpara" in Diffusion Map is not suitable(Usually too large). Please Change. ') 
%     % end
%     if (max(K(:,1))-min(K(:,1)))<0.7    
%         warning('Kernel Matrix varies very limit. Parameter "options.kpara" in Diffusion Map is not suitable(Usually too large). Please Change. ') 
%     end
% 
%     %------------------------------------------
%     if R~=R2 
%         warning('R~=R2') 
%     end
% 
%     R2*LamdaS*L'-P;
% %------------------------------


%--------------------------------------------------------------------------
if ~isreal(LamdaS(1:dim_new+1,1:dim_new+1))
     manifoldLearning('Eigenvalue (of P) used involves complex number') 
end
if ~isreal(R(:,1:dim_new+1))
     manifoldLearning('Eigenvector (of P) used involves complex number') 
end

% END of Check point!------------------------------------------------------

% ------------------------------------------
% % Truncation/Dimension Reduction
% V=V(:,1:Dim_new);
% Lamda=Lamda(1:Dim_new,1:Dim_new);
% 
% Vn=Vn(:,1:Dim_new);
% Lamdan=Lamdan(1:Dim_new,1:Dim_new);

%------------------------------------------
% Diffusion Process and new coordinate systemn

lamda_t=LamdaS.^options.t;                    %Diffusion process by step t
% Z=lamda_t(2:dim_new+1)*(V(:,2:dim_new+1))';

% Z=R(:,2:dim_new+1)*lamda_t(2:dim_new+1,2:dim_new+1);
% The information in first direction is void as it is a constant. Do not
% use. Do it in this way for easy use with preimage solution.
Z=R(:,1:dim_new+1)*lamda_t(1:dim_new+1,1:dim_new+1);

% Z2=R(:,2:dim_new+1);      % Another way of giving result


switch options.Ztype
    case 0  %without 1st component
        model.Z1=Z(:,1);    
        model.Z11=Z(1,1);
        Z=Z(:,2:end);          
    case 1  %with 1st component
        Z=Z;
        model.Z11=Z(1,1);
end
        
    
%% Recording & Output
%-----------------------------------------
% Basic Parameter Information (Input Information)

% model.metric=options.metric;
% model.kernel=options.kernel;
% model.kpara=options.kpara;
% model.dim_new=dim_new;
% model.t=t;
model.Z=Z;


model.options=options;
model.X=X;

model.eigenvalues=diag(LamdaS);

model.R=R(:,1:dim_new+1);                        %  Dimension Reduced right Eigenvectors of P
model.L=L(:,1:dim_new+1);                        %  Dimension Reduced left  Eigenvectors of P
model.S=S(:,1:dim_new+1);                        %  Dimension Reduced Eigenvectors of P'
model.Lamda=LamdaS(1:dim_new+1,1:dim_new+1);      %  Dimension Reduced Eigenvalue of P and P' (They are the same)

%-----------------------------------------
% Full Information. Conditional output. (For Further Research without recalculation OR Error Report. Mass Memoey needed)
if options.FullRec == 1
    
    %-----------------------------------------
    % Processed Information (Result Information)
    model.K=K;
    model.D=D;              % Degree Matrix of K

    model.P=P;
    model.Pprime=Pprime;              % Degree Matrix of K
    model.Lt=P^options.t;
    model.dist=Distance;
    
    model.R_full=R;              % Full Eigenvectors of L  (The Laplacian Matrix)
    model.L_full=L;            % Full Eigenvectors of Ln (The Normalized Laplacian Matrix)
    model.S_full=S;
%     model.Lamd_full=Lamda;      % Full Eigenvalue of L and Ln (They are the same)
%     model.LamdaS_full=LamdaS;      % Full Eigenvalue of L and Ln (They are the same)

end

return


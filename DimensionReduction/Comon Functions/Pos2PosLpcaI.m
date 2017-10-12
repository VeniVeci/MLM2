function [X_star]= Pos2PosLpcaI(X,Z,Dist,Z_star,options)
% Positions (in feature space) to positions (in original space)
%
% Synopsis:
% [X_star]= Pos2PosLpcaI(X,Z,Dist,Z_star,options)
%
% Description:
% Positions in feature space to positions coordinate calculation using Different
% Schemes 
% 
% Input:
%        X                         % [sample x dimension] observation in original space
%        Z                         % [sample x dim_fea] corresponding feature observations
%        dist                      % [num_Z_star x sample] Distance between X and
%                                    X_star in original space.
%        Z_star                    % [num_Z_star x dim_fea] feature observations to be
%        options                   % parameters         
%
% Output:
%       X_star                    % [num x dimension] recovered observation in original space 
% 
%       model [structure]         % Information about the original data and model
%            .options             % Function parameters
%
%
% Example:
% See also:
% About: 
% Pcakage Require:
% See also:
% 
% Modifications:
% WeiX, 1-11-2016, DayI Create 
%
%% Initialization and Parameters
[num_X,DimOrig]=size(X);
[num_Z_star,Dim_Z_star]=size(Z_star);     %Dimension of feature space
[num_Z,Dim_Z]=size(Z); 

if Dim_Z ~= Dim_Z_star      %Identify Dimension of feature space
    warning('Dimension dismatch. Auto Truncation Actived')
    Z=Z(:,1:Dim_Z_star);
    DimF=Dim_Z_star;
else 
    DimF=Dim_Z_star;
end

%Check for input distances
if (min(Dist(:))<0)
    error('Distances must be positive value');
end
    
% if nargin < 2, options = []; else options=c2s(options); end
if nargin < 5, options = []; end
if ~isfield(options,'dim_new'), options.dim_new = DimF; end                     % Default output information;
if ~isfield(options,'InMethod'), options.InMethod = 'LSE'; end              % Default Interpretation method;
if ~isfield(options,'FullRec'), options.FullRec = 0; end                    % Default output information;
if ~isfield(options,'neighbor'), options.neighbor = num_X; end              % Default kernel function

% %Rearrange order
% if options.neighbor ~= num_X
%     [Dist,index]=sort(Dist);
%     Dist = Dist(1:options.neighbor,:);
%     X = X(index(1:options.neighbor),:);
% end


%% Main 
for i=1:num_Z_star
    
    [Disti,index]=sort(Dist(i,:));
    Xi = X(index(1:options.neighbor),:);
    Zi = Z(index(1:options.neighbor),:);
    
    %PCA process via SVD on X
    dim_new=options.dim_new;
%     [U,S,V]=svd(Xi);                     % X=U*S*V'
    [U,S,V]=svd(Xi,'econ');                % quick SVD. X=U*S*V'
    
    [num_Xi,dim_Xi]=size(Xi);
    if (dim_new>num_Xi)|(dim_new>dim_Xi)
        error('Not enough neighbors for svd-econ');
    end
    
    U_temp=U(:,1:dim_new);
    S_temp=S(1:dim_new,1:dim_new);      
    V_temp=V(:,1:dim_new);              % X~=U_temp*S_temp*V_temp'

    Yi=U_temp*S_temp;                    % Projections
    
    %% Interprete from Z* to Y*
    switch options.InMethod
        case 'LSE'            
            A=Zi\Yi;              %Linear approximation Z to Y via LSE Y=ZA
            Y_stari=Z_star(i,:)*A;
%             X_star(i,:)=Y_stari*V_temp';
                 
        case 'ANN'
            %% Multi-variate ANN
            net=feedforwardnet(10); % One hidden layer with nn nodes; for more layers, 
            % use [nn1 nn2 nn3 ... nnJ] for J layers with nnj nodes in the jth layer 
            net = init(net); % Reinitialize weights
            net.divideParam.trainRatio=0.9; % Fraction of data used for training (cross-validation)
            net.divideParam.valRatio=(1-net.divideParam.trainRatio);% /2; % Fraction of data used for validation
            net.divideParam.testRatio=0;% (1-net.divideParam.trainRatio)/2; % Fraction of data used for testing
            [net,tr] = trainlm(net,Zi',Yi'); % Feedforward with Levenberg-Marquardt backpropagation
%             [net,tr] = trainbr(net,Zi',Yi'); % Bayesian Regulization
            Y_stari=net(Z_star(i,:)');   
            Y_stari=Y_stari';
%             X_star(i,:)=Y_stari*V_temp';
                
        otherwise 
            error('No such Interpret Method')
    end
    
    X_star(i,:)=Y_stari*V_temp';
    
    
end





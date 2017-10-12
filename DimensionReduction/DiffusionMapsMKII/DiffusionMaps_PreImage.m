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
%  WeiX, 4-1-2016,  Update code for updated DiffusionMap with Ztype.
%% ----------------------Initialization------------------------------------
if nargin <= 2, options = []; end


[Num_Zstar,Dim_Zstar]=size(Z_star);
[Num_X,Dim_X]=size(model.X);

X_star=zeros(Num_Zstar,Dim_X);          % Assign memory
Dist_star=zeros(Num_Zstar,Num_X);       % Assign memory

Lamda=model.Lamda;
L=model.L;
R=model.R;

%Dimension dismatch
if Dim_Zstar~=model.options.dim_new+1
    warning('Dimension dismatch. Program automatically truncate to match');
end

%Find k_**
switch model.options.kernel
    case 'gaussian'
        k_starstar = exp(-(0).^2/model.options.kpara);    % Actually the "0" could change due to differnet metric method.       
    otherwise
        error('Error: Undefined type of kernel function.');    
end

%Check Ztype
switch model.options.Ztype
    case 0  %without 1st component
        Z_star=[repmat(model.Z11,[Num_Zstar,1]),Z_star];
    case 1  %with 1st component
        Z_star1=mean(Z_star(:,1));
        if (Z_star-model.Z11)>1e-10
            warning('Z_star 1st component do no match with model')
        end
end        
[Num_Zstar,Dim_Zstar]=size(Z_star); %Update Z_star information


%% -----------------------Main--------------------------------
% R(:,1)*lamda_t(1,1)

for i = 1:Num_Zstar

%     r(i,:)=Z_star(i,:)*(Lamda^(-model.options.t));    
%     P_stardot=r(i,:)*Lamda*L';
%     P_stardot2=Z_star(i,:)*L';

% Use to deal with Dim_Zstar~=model.options.dim_new+1
    r(i,:)=Z_star(i,:)*(Lamda(1:Dim_Zstar,1:Dim_Zstar)^(-model.options.t));    
    P_stardot=r(i,:)*Lamda(1:Dim_Zstar,1:Dim_Zstar)*L(:,1:Dim_Zstar)';
         
    
    A=P_stardot'*ones(1,Num_X)-eye(Num_X);
    b=-k_starstar*P_stardot';
%     b=k_starstar*P_stardot';
%     k_star=A\b;
    
    lb=zeros(Num_X,1);
    ub=ones(Num_X,1);
%     k_star = lsqlin(A,b,[],[],[],[],1e-5,1);
    k_star = lsqlin(A,b,[],[],[],[],lb,ub);
    %--------------
%     k_starorig=model.K(i,:)';
%     nums=200-1
%     [U,S,V]=svd(A);
%     k_star2=V(:,1:nums)*S(1:nums,1:nums)*U(:,1:nums)'*b;
%     
%     
%     % Some Fix 
%     k_star=abs(k_star);
%     
%     sumkstar=sum(k_star)
%     power=round(log10(sumkstar));
%     k_star=k_star./10^power
%     sumkstar=sum(k_star)+k_starstar
%     
%     k_starorig=model.K(i,:)';
%     p_starrecover=k_star./sumkstar; % used to compare P_stardot
%     
%     k_compare=k_star./k_starorig;
%     

    %Monitor Point!
    K_rec(i,:)=k_star;
           
%     k_star = lsqlin(A,b,[],[],[],[],0,1);
    
%     plot([p_starrecover';P_stardot]');
    
    % Test 
    % real P_stardot
%     P_stardotorig=model.P(i,:)
%     
%     A=P_stardotorig'*ones(1,Num_X)-eye(Num_X);
%     A=P_stardot'*ones(1,Num_X);
%     b=-k_starstar*P_stardotorig';
% %     b=k_starstar*P_stardot';
%     k_starorig=A\b;
    

%     Z_star_complete=[v11*Lamda(1,1),Z_star(i,:)];  
% %     v_star=Z_star(i,:)*Lamda(1:end,1:end)^(-model.options.t);
%     v_star=Z_star_complete*Lamda(1:end,1:end)^(-model.options.t);
%     
%     v_starM1=Z_star(i,:)*Lamda(2:end,2:end)^(-model.options.t);  %V_starM1= Vectorc_star minus the ist component
%     v_star=[v11,v_starM1];
%     
%     d_star=sqrt(k_starstar)/sqrt(v_star*C^2*Lamda*v_star'); 
%     % k_starstar=d_star*v_star*C^2*Lamda*v_star'*d_star;
%     k_star=D*V*C^2*Lamda*v_star'*d_star;
%     
%     K_orig=D*V*C^2*Lamda*V'*D;

% -------------!!!!!!!!! Correction !!!!!!!!!!!!-------------------------- 
% k_star=(k_star>=0).*k_star;             % Ensure k_star is positive.
% k_star=k_star+(k_star==0).*1e-50;       % Ensure 0 in k_star are limt to 0 Rather than real '0'..

    if max(k_star)>1
        warning('k_star contains value larger than 1. Force auto correction')
        k_star=k_star./max(k_star);              % Ensure maximum value is 1

        %Report Monitor Point! 
        K_rec(i,:)=k_star;  
    end
% ------------------------------------------------------------------------

    switch model.options.kernel
        case 'gaussian' 
            dist_star= sqrt(-log(k_star)*model.options.kpara);
        otherwise
            error('Error: Undefined type of kernel function.');    
    end
    
%OLD VERSION    
%     X_stari= Dist2pos(model.X,dist_star,options);
%     X_star(i,:)= X_stari;
    
    switch options.type
        case 'LpcaI'
%             options.InMethod = 'LSE';          		 % Default Interpretation method;            
%             X_star(i,:)= Pos2PosLpcaI(model.X,model.Z,dist_star',Z_star(i,2:end),options); %%%%%%%+++ERROR
            switch model.options.Ztype
                case 0
                    X_star(i,:)= Pos2PosLpcaI(model.X,model.Z,dist_star',Z_star(i,2:end),options); 
                case 1
                    X_star(i,:)= Pos2PosLpcaI(model.X,model.Z,dist_star',Z_star(i,:),options); 
            end      
            
            
        otherwise
            X_star(i,:)= Dist2pos(model.X,dist_star,options);
    end
    
    
    
    
end
    


return % End of Function



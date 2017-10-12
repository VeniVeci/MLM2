%Temp Fix problem form P to K in diffusion maps

clear 

Num_data=100;

%% Data Generation Swiss Roll Dataset
t1=rand(Num_data,1)*4*pi;   % Theater
t2=rand(Num_data,1)*20;     % change 20 to 10 will show different
t1=sort(t1);                
X(:,1)=t1.*cos(t1);         % X
X(:,2)=t2;                  % Y
X(:,3)= t1.*sin(t1);        % Z
Xcolor = t1;                     % Color according to Theater
Xsize = ones(Num_data,1)*10;    % Size. Constant





% DiffusionMap options
options.metric ='euclidean';
options.kernel ='gaussian'; 
options.kpara = 10000;    
options.kAuto=1;

          
options.t = 1;                     
options.FullRec = 0;      

options.dim_new = 3;    
dim_new=options.dim_new;

%% Diffusion maps 
[num,dim]=size(X);

Distance =pdist2(X,X,options.metric);

if options.kAuto==1
    options.kpara=sum(Distance(:).^2)/(num^2);
end

% Calculating the Kernel Matrix
switch options.kernel
    case 'gaussian'
        K = exp(-Distance.^2/options.kpara);            
    otherwise
        error('Error: Undefined type of kernel function.');    
end

d=sum(K,2);                             % row sum
D=diag(d);                              % Degree Matrix. It is a diagonal matrix. 
P=D^(-1)*K;                             % The Markov matrix
Pprime=D^(-0.5)*K*D^(-0.5);             % The Normalized Markov matrix

%-----------------------------------------
% % Eigen Decomposition of The Markov matrix
[R2,Lamda] = eig(P);  
[Lamda,Index] = sort(diag(Lamda),'descend'); % sort eigenvalues in descending order
R2 = R2(:, Index);                             % sort eigenvector in corresponding order to its eigenvalue
Lamda=diag(Lamda);

%------------------------------------------
% % Eigen Decomposition of The Normalized  Markov matrix
[S,LamdaS] = eig(Pprime);  
[LamdaS,Index] = sort(diag(LamdaS),'descend');
S = S(:, Index);
LamdaS=diag(LamdaS);                     % LamdaS must = Lamda

%------------------------------------------
% % Result
R=D^(-0.5)*S;  % R2 must =R     % Right eigenvalue of P
% 
% if R~=R2 
%     warning('R~=R2') 
% end
% 
L=D^(0.5)*S;
% 
% R2*LamdaS*L'-P;
% 
lamda_t=Lamda.^options.t;                    %Diffusion process by step t

% Z=R(:,2:dim_new+1)*lamda_t(2:dim_new+1,2:dim_new+1);
% The information in first direction is void as it is a constant. Do not
% use. Do it in this way for easy use with preimage solution.


Z=R(:,1:dim_new+1)*lamda_t(1:dim_new+1,1:dim_new+1);


X2=[X;X(1,:)];

Distance =pdist2(X2,X2,options.metric);

% if options.kAuto==1
%     options.kpara=sum(Distance(:).^2)/(num^2);
% end

% Calculating the Kernel Matrix
switch options.kernel
    case 'gaussian'
        K2 = exp(-Distance.^2/options.kpara);            
    otherwise
        error('Error: Undefined type of kernel function.');    
end

d2=sum(K2,2);                             % row sum
D2=diag(d2);                              % Degree Matrix. It is a diagonal matrix. 
P2=D2^(-1)*K2;                             % The Markov matrix
Pprime2=D2^(-0.5)*K2*D2^(-0.5);             % The Normalized Markov matrix




%-----------------------------------------------------------------------
% Test Zone

%     P_stardot=r(i,:)*Lamda*L';
%     P_stardot2=Z_star(i,:)*L';

    P_stardot_orig=P2(end,1:end-1);
%     P_stardot=P_stardot(:,1:end-1);
%     P_stardot=P_stardot+rand(1,num)*0.05;
    P_stardot_approx=Z(1,:)*L(:,1:dim_new+1)';
    
    % Here shows that P_stardot_approx is simliar to P_stardot_orig
    P_compare=P_stardot_orig./P_stardot_approx;
    
    
    P_stardot=P_stardot_orig;
    P_stardot=P_stardot_approx;
    
    A=P_stardot'*ones(1,num)-eye(num);
%     b=-k_starstar*P_stardot';
    b=-1*P_stardot';

%     k_star=A\b;
    k_star = lsqlin(A,b,[],[],[],[],0.000001,1);
    
    %--------------
    k_starorig=K2(end,1:end-1)';
    
    %Here shows that, however, k_star and k_starorig are completely
    %different
    k_compare=k_starorig./k_star;
%     k_compare=k_compare*10e15;
    
    
    nums=200-1
    [U,S,V]=svd(A);
    k_star2=V(:,1:nums)*S(1:nums,1:nums)*U(:,1:nums)'*b;
        
    % Some Fix 
    k_star=abs(k_star);
    
    sumkstar=sum(k_star)
    power=round(log10(sumkstar));
    k_star=k_star./10^power
    sumkstar=sum(k_star)+k_starstar
    
    k_starorig=model.K(i,:)';
    p_starrecover=k_star./sumkstar; % used to compare P_stardot
    
    k_compare=k_star./k_starorig;

%-----------------------------------------------------------------------






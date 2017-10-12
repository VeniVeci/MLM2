function [Y_star]= ANN_prediction_crossvalidation(X,Y,X_star,Kfold)
% Artificial neural network prediction with crossvalidation
%  
% Synopsis:
%  [Y_star]= ANN_prediction_crossvalidation(X,Y,X_star,Kfold)
%
% Description:
% predict the y_star distribution using given information through ANN.
% !!There are some parameter one can change in the function.!!!
% 
% steps for the program:
% 
% Input:
%  X        [dimension_X x samples] indicates Training data X.
%  Y        [dimension_X x samples] indicates Training data Y.
%  X_star   [dimension_Y x test points] indicates conditions of prediction Y_star
%  Kfold    [1 x 1] value unmber indicates fold used for crossvalidation
% 
% Output:
% Y_star      [dimension_Y x test points] mean values of predictions
% Y_var_star  [dimension_Y x test points] variance of preditions
%
% Example:
%
% See also 
% 
% About: 
% 
% Modifications:
% 19-jul-2013, WeiX, first edition 


[Dim_X,Num_X]= size(X);
[Dim_Y,Num_Y]= size(Y);
[Dim_Xstar,Num_Xstar]= size(X_star);


if Kfold>1
    
    Num_inFold=fix(Num_X/Kfold);
    SSE_min=inf; %can be adjusted according to dataset 
    
    for fold=1:Kfold

        X_star_temp=X(:,(fold-1)*Num_inFold+1:fold*Num_inFold);
        X_temp= (removerows(X','ind',(fold-1)*Num_inFold+1:fold*Num_inFold))'; 

        Y_starorig_temp=Y(:,(fold-1)*Num_inFold+1:fold*Num_inFold);
        Y_temp= (removerows(Y','ind',(fold-1)*Num_inFold+1:fold*Num_inFold))'; 

        [Dim_X_temp,Num_X_temp]= size(X_temp);
        [Dim_Y_temp,Num_Y_temp]= size(Y_temp);
        [Dim_Xstar_temp,Num_Xstar_temp]= size(X_star_temp);


        net=feedforwardnet(10); % One hidden layer with nn nodes; for more layers, 
        %net=cascadeforwardnet(10);
        % use [nn1 nn2 nn3 ... nnJ] for J layers with nnj nodes in the jth layer 
        net = init(net); % Reinitialize weights
        net.divideParam.trainRatio=0.9; % Fraction of data used for training (cross-validation)
        net.divideParam.valRatio=(1-net.divideParam.trainRatio);% /2; % Fraction of data used for validation
        net.divideParam.testRatio=0;% (1-net.divideParam.trainRatio)/2; % Fraction of data used for testing
        
        %Stop nntriantool popping up
        net.trainParam.showWindow=0; %default is 1
        
        
        
        [net,tr] = trainlm(net,X_temp,Y_temp); % Feedforward with Levenberg-Marquardt backpropagation    

        for i=1:Num_Xstar_temp        
        Y_star_temp(:,i)=net(X_star_temp(:,i));   
        end

        SquErr=(Y_starorig_temp-Y_star_temp).^2;
        SSE=sum(SquErr(:));

        if SSE<SSE_min
            SSE_min=SSE;
            net_best=net;
        end
        
    end
    
    for i=1:Num_Xstar
        Y_star(:,i)=net_best(X_star(:,i));    
    end   

    return    
    
    
    
else % when there is no cross-validation (Kfold=1)
   
    net=feedforwardnet(10); % One hidden layer with nn nodes; for more layers, 
    % use [nn1 nn2 nn3 ... nnJ] for J layers with nnj nodes in the jth layer 
    net = init(net); % Reinitialize weights
    net.divideParam.trainRatio=0.9; % Fraction of data used for training (cross-validation)
    net.divideParam.valRatio=(1-net.divideParam.trainRatio);% /2; % Fraction of data used for validation
    net.divideParam.testRatio=0;% (1-net.divideParam.trainRatio)/2; % Fraction of data used for testing
    [net,tr] = trainlm(net,X,Y); % Feedforward with Levenberg-Marquardt backpropagation


            %current1=net(xtest);
            %current2=net(xtest2);

    for i=1:Num_Xstar        
        Y_star(:,i)=net(X_star(:,i));    
    end
    
    return
           
end


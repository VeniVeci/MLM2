function [Y_star]= ANN_prediction(X,Y,X_star)
% Artificial neural network prediction
%  
% Synopsis:
%  [Y_star]= ANN_prediction(X,Y,X_star)
%
% Description:
% predict the y_star distribution using given information through ANN.
% There are some parameter one can change in the function.!!!
% 
% steps for the program:
% 
% Input:
%  X        [dimension_X x samples] indicates Training data X.
%  Y        [dimension_X x samples] indicates Training data Y.
%  X_star   [dimension_Y x test points] indicates conditions of prediction Y_star
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


% X=X';
% Y=Y';
% X_star=X_star';

[Dim_X,Num_X]= size(X);
[Dim_Y,Num_Y]= size(Y);
[Dim_Xstar,Num_Xstar]= size(X_star);


net=feedforwardnet(10); % One hidden layer with nn nodes; for more layers, 
% use [nn1 nn2 nn3 ... nnJ] for J layers with nnj nodes in the jth layer 
net = init(net); % Reinitialize weights
net.divideParam.trainRatio=0.9; % Fraction of data used for training (cross-validation)
net.divideParam.valRatio=(1-net.divideParam.trainRatio);% /2; % Fraction of data used for validation
net.divideParam.testRatio=0;% (1-net.divideParam.trainRatio)/2; % Fraction of data used for testing
% [net,tr] = trainlm(net,X,Y); % Feedforward with Levenberg-Marquardt backpropagation
[net,tr] = trainbr(net,X,Y); % Feedforward with Levenberg-Marquardt backpropagation
        
        %current1=net(xtest);
        %current2=net(xtest2);
        
for i=1:Num_Xstar        
    Y_star(:,i)=net(X_star(:,i));    
end
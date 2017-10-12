%Test_PreImage
%Comparing reconstruction performance between different metohds and
%parameters
clear
%-----------------------------------------------------------------------

%Assign parameter
new_dim_start=7;    
new_dim_end  =7;
kfold=10;                       %Defining cross-validation
iteration=10;                    %Defining test times

num_train=100;                  %Defining the input
num_test=1;

%-----------------------------------------------------------------------
%temp variable
Index_Fig=1;
Rec=zeros(1,6);

%-----------------------------------------------------------------------
%Getting input

    [X,Y,X_star,Y_starorig]=DataGenerate(num_train,num_test);

        %X=load('Inputs2.txt');
        %Y=load('Outputs2.txt');
        %X_star=load('InputTest14');
        %Y_starorig=load('OutputTest14');

    %---------------------------------------------------------------
    [np_train,Dim_X]=size(X);
    [np_test, Dim_Y]=size(Y_starorig);

    X_whole=[X;X_star];
    Y_whole=[Y;Y_starorig];

    [X_whole,X_digi]=DataDigiRescale(X_whole);
    X=X_whole(1:np_train,:);
    X_star=X_whole(np_train+1:end,:);

    [Y_whole,Y_digi]=DataDigiRescale(Y_whole);
    Y=Y_whole(1:np_train,:);
    Y_starorig=Y_whole(np_train+1:end,:);

    %--------------------------------------------------------------

    X=X';
    Y=Y';
    X_star=X_star';
    Y_starorig=Y_starorig';

    %------------------------------------------------------------------------------------------
    % Testing
    j=1;
    for i = -7:1:-7
    options = struct('ker','sigmoid','arg',[10^i,0],'new_dim',7); %0.000001
    Rec_SSE(:,j)=PreImageTest_KPCA(Y,options,new_dim_end);
    j=j+1;
    end

    
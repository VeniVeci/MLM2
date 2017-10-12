% Sobo_Generator_Demo
%
%  Modification:
%  WeiX, Nov 26th 2014, First Edition

Num=1000;
rang=[5,50;-5,-50];
    
[X] = Sobo_Generator(Num,rang);

figure(1)
scatter(X(:,1),X(:,2));

Z = exp ( (-X(:,1).^2 - X(:,2).^2)/100 );
figure(2)
scatter3(X(:,1),X(:,2),Z)
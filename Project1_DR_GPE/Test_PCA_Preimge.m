% PCA preimge and Principal direction (vector) reconstruction using X and Z
% Experimental project.

% Require Dr PCA package


clear 

x=rand(100,20);
[z,model] = Lpca(x,3);
 z2=x*model.Key';

for i=1:100
    
    n=i;
    xx=x(1:n,:);
    zz=z(1:n,:);
%      kk=zz\xx;
     kk=xx\zz;
    
%     reckk(i,:)=kk;
end



    n=20;
    xx=x(1:n,:);
    zz=z(1:n,1);

    kk=zz'*inv(xx);
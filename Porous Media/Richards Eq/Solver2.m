%% Richards equation Demo script 

nNode=10;
nTime=100;.

deltaZ=0.1;
deltaTime=1;

nMaxIteration=100;

timeStep=0.01;

H=zeros(nNode,nTime);

H_temp=zeros(nNode,nMaxIteration);

%%% 
H=zeros(nTime,nNode,nMaxIteration);
K=zeros(nTime,nNode,nMaxIteration);

time=[0:deltaTime:100];










for t=1:length(time)
    
    
    for j=2:nMaxIteration
        
        for i = 2:nNode-1
    
            K=1;
            a=(K(i)+K(i-1))/(-2*deltaZ^2);
            b=c(i)/deltaTime+(K(i+1)+2*K(i)+K(i-1))/(2*deltaZ^2);
            d=-(K(i+1)+K(i))/(2*deltaZ^2);
            e=(-K(i+1)+K(i-1))/(2*deltaZ)+c(i)*H(i,t)/deltaTime;

        %     A(i-1,i-2)=a;
        %     A(i-1,i-1)=b;
        %     A(i-1,i)=d;

            A(i,i-1)=a;
            A(i,i)=b;
            A(i,i+1)=d;

            b(i)=e;

        end
        h=b\A;
        
        
        
        
        
        
        
        a=K(t,i,j)
        
        a=(K(i)+K(i-1))/(-2*deltaZ^2);
        
        
        
        
        
    end
    
    
   
    
    
    
    
    
end



t=1;
for i = 2:nNode-1
    
    left=
    a=(k(i)+k(i-1))/(-2*deltaZ^2);
    b=c(i)/deltaTime+(k(i+1)+2*k(i)+k(i-1))/(2*deltaZ^2);
    d=-(k(i+1)+k(i))/(2*deltaZ^2);
    e=(-k(i+1)+k(i-1))/(2*deltaZ)+c(i)*H(i,t)/deltaTime;
    
%     A(i-1,i-2)=a;
%     A(i-1,i-1)=b;
%     A(i-1,i)=d;
    
    A(i,i-1)=a;
    A(i,i)=b;
    A(i,i+1)=d;
    
    b(i)=e;
   
end
h=b\A;



function [h]= picardIteration(h0,nNode,deltaZ,deltaTime)

k=1;
h=h0;

%CSC approximation 





for j=1:nMaxIteration    
    
    for i = 2:nNode-1
        
        
        %theata
        theata=alpha*(theata_s-theata_r)/(alpha+abs(h)^beta+theata_r);
        
        %k
        k=k_s*rho/(rho+abs(h)^r);
        
        %C
        c(i)=(theata(i)-theata0(i))/(h(i)-h0(i));
        
        
        
        a=(k(i)+k(i-1))/(-2*deltaZ^2);
        b=c(i)/deltaTime+(k(i+1)+2*k(i)+k(i-1))/(2*deltaZ^2);
        d=-(k(i+1)+k(i))/(2*deltaZ^2);
        e=(-k(i+1)+k(i-1))/(2*deltaZ)+c(i)*H(i,t)/deltaTime;

    %     A(i-1,i-2)=a;
    %     A(i-1,i-1)=b;
    %     A(i-1,i)=d;

        A(i,i-1)=a;
        A(i,i)=b;
        A(i,i+1)=d;

        b(i)=e;

    end
    h=b\A;
    
end

    



end 

function result=theata(h)
theata_s=1;
theata_r=2;
alpha=1;
beta=2;
theata_r=1;

result=alpha*(theata_s-theata_r)/(alpha+abs(h)^beta+theata_r);

end






function X=LESsolver_Compare(A,b)


%     option(1)=  'NGE ';
%     option(2)=  'TDMA';
%     option(3)=  'BSOR';
%     option(4)=  'GE  ';
%     option(5)=  'JCB ';
%     option(6)=  'GS  ';
%     option(7)=  'SOR ';
%     option(8)=  'CG  ';



    S = ['NGE ';'TDMA';'BSOR';'GE  ';'JCB ';'GS  ';'SOR ';'CG  '];
    S = ['NGE ';'TDMA';'BSOR';'GE  ';'CG  ']
    

    X(:,1)=A\b;
    for i=1:5
             
        option = strtrim(S(i,:));
        X(:,i+1)=LESsolver(A,b,option);
        
    end
    
    
end

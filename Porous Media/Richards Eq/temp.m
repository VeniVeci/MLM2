H(:,1)=h_init;
h=h_init;
for t=1:nTime
    
    h_PreviousTime=h;
    h0=h_PreviousTime; %???
    
    for j=1:nMaxIteration     

        %Assemble matrix for A_l*x=b_l
        matrixA=zeros(nNode);
        columnB=zeros(nNode,1);
        
        for i = 2:nNode-1
                
            %c by SCS approximation
%             h0(i)=h0(i)+0.0000001;    %To avoid NaN      
            c(i)=(theata(h0(i))-theata(h_PreviousTime(i)))/(h0(i)+0.00000001-h_PreviousTime(i)); % add 0.0000001 To avoid NaN
            
            %c by analytical derivatives
%             c(i)=theataDif(h(i));
            
            
            a=(k(h0(i))+k(h0(i-1)))/(-2*deltaZ^2);
            b=c(i)/deltaTime+(k(h0(i+1))+2*k(h0(i))+k(h0(i-1)))/(2*deltaZ^2);
            d=-(k(h0(i+1))+k(h0(i)))/(2*deltaZ^2);
            e=(-k(h0(i+1))+k(h0(i-1)))/(2*deltaZ)+c(i)*h_PreviousTime(i)/deltaTime;

            matrixA(i,i-1)=a;
            matrixA(i,i)=b;
            matrixA(i,i+1)=d;

            columnB(i)=e;

        end
        
        %Deal with B.C.(truncate the system matrix and take care of the
        %consequence)
        
        columnB=columnB-h_init(1).*matrixA(:,1);
        columnB=columnB-h_init(end).*matrixA(:,end);
        
        
        matrixA=matrixA(2:nNode-1,2:nNode-1);
        columnB=columnB(2:nNode-1,1);
        
%         columnB(1)=columnB(1)-matrixA(1,2)*h_init(1);
%         columnB(end)=columnB(end)-matrixA(end-1,end)*h_init(end);
        
       
        %Utilize sparsity
%         sparseA=sparse(matrixA);
%         h=sparseA\columnB;
        
        
        %Solve linear system
        h=matrixA\columnB;
        
        %Complete h with B.C.
        h=[h_init(1);h;h_init(end)];
        
        
        smeIte=sum((h-h0).^2);
        if sqrt(smeIte)<miniIteError 
            break 
        else
            h0=h;
        end
        

    end
    
    H(:,t+1)=h;
    
end
computeTime=toc;






%% SAVE FOR TEMP
%             %BC point
%             %top
%             if k==1 
%             %
%             weightTop=0;
%             K_top=0;
%                 
%             %bottom
%             if k==nNodeZ
%                 
%                 weightBottom=0;              
%                 K_bottom=0;       
%                 
%             %left
%             if j==1 
%                 weightLeft=0;
%                 K_top=0; %if k=1
% 
%             %right
%             if j==nNodeX
%             
%             %inner
%             else 
%                 weightLeft=-deltaTime/deltaX^2*D_left;
%                 weightRight=-deltaTime/deltaX^2*D_right;
%                 
%                 weightTop=-deltaTime/deltaZ^2*D_top;
%                 weightBottom=-deltaTime/deltaZ^2*D_bottom;
%                 
%                 weightCenter=1 + deltaTime/deltaX^2*D_left ...
%                                + deltaTime/deltaX^2*D_right ...
%                                + deltaTime/deltaZ^2*D_top ...
%                                + deltaTime/deltaZ^2*D_bottom;
%                            
%                 b=theata0-deltaTime/deltaZ*(K_bottom-K_top);
%             end

%%



%    for j=1:nNodeX
%         for k=1:nNodeZ
%             
%             index=index+1; %next grid point 
%             
% %             %Coordinate to index
% %             indexNextX=index+1;
% %             indexLastX=index-1;
% %             
% %             indexNextZ=index+nNodeX+1;
% %             indexLastZ=index-nNodeX-1;
%             
%             
%             %--
%             %check if top BC
%             if k==1
%                 D_HalfLastZ=0;
%                 K_HalfLastZ=0;
%                 indexLastZ =0;
%             else 
%                 D_HalfLastZ=(D(k,j)+D(k-1,j))/2; 
%                 K_HalfLastZ=(K(k,j)+K(k-1,j))/2;
%                 indexLastZ =index-1;
%             end
% 
%             %check if bottom BC
%             if k==nNodeZ
%                 D_HalfNextZ=0;
%                 K_HalfNextZ=0;
%                 indexNextZ =0;
%             else 
%                 D_HalfNextZ=(D(k,j)+D(k+1,j))/2; 
%                 K_HalfNextZ=(K(k,j)+K(k+1,j))/2;
%                 indexNextZ =index+1;
%             end
% 
%             %check if left BC
%             if j==1
%                 D_HalfLastX=0;
%                 indexLastX =0;
%             else 
%                 D_HalfLastX=(D(k,j)+D(k,j-1))/2; 
%                 indexLastX =index-nNodeZ;
%             end
% 
%             %check if right BC
%             if j==nNodeX
%                 D_HalfNextX=0;
%                 indexNextX=0;
%             else 
%                 D_HalfNextX=(D(k,j)+D(k,j+1))/2; 
%                 indexNextX =index+nNodeZ;
%             end       
%             
%             
%             %Assemble system 
%             weightNextX=-deltaTime/deltaX^2*D_HalfNextX;
%             weightLastX=-deltaTime/deltaX^2*D_HalfLastX;
%             
%             weightNestZ=-deltaTime/deltaZ^2*D_HalfNextZ;
%             weightLastZ=-deltaTime/deltaZ^2*D_HalfLastZ;
%             
%             weightSelf = 1 - weightNextX ...
%                            - weightLastX ...
%                            - weightNestZ ...
%                            - weightLastZ ;
%             
%                       
%             b(index)=H_PreviousTime(k,j)-deltaTime/deltaZ*(K_HalfNextZ-K_HalfLastZ);
%             
%             if indexNextX>0 A(index,indexNextX)=weightNextX; end
%             if indexLastX>0 A(index,indexLastX)=weightLastX; end
%             if indexNextZ>0 A(index,indexNextZ)=weightNestZ; end
%             if indexLastZ>0 A(index,indexLastZ)=weightLastZ; end
% 
%         end
% 
%     end


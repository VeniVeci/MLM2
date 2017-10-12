
clear 


h1 = waitbar(0,'Loop1');


for i = 1:10
    
%     close(h1)  
%     h1 = waitbar(0,'Loop1');
%     waitbar(i/10,h1,'Loop1')

    waitbar(i/10,h1)
    
    h2 = waitbar(0,'Loop2');
    for j=1:10
        
        waitbar(j/10,h2)
    
    end
    close (h2)   

      
      
      
end
      
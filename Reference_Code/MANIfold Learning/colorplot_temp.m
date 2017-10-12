        N=1000;

        
        tt = (3*pi/2)*(1+2*rand(1,N));  
        height = 21*rand(1,N);
        handles.X = [tt.*cos(tt); height; tt.*sin(tt)]';
        handles.ColorVector = tt';

        figure(2);
        scatter3(handles.X(:,1),handles.X(:,2),handles.X(:,3),12,handles.ColorVector,'filled');
        
        figure(3);
        scatter3(handles.X(:,1),handles.X(:,2),handles.X(:,3),handles.ColorVector,handles.ColorVector,'filled');
function [Data] = Data_GeneratorM(N,type,options)
%%  Data_Generator M (mesh)
%   dataset generation for test or further use. Data uniformly distributed.
%
% Input: 
%       N              --- Number of data point
%       type           --- Type of dataset. 
%                             1. 'SwissRoll'
%                             2. 'SwissHole'
%                             3. 'CornerPlanes'
%                             4. 'PuncturedSphere'
%                             5. 'TwinPeaks'
%                             6. '3DClusters'
%                             7. 'ToroidalHelix'
%                             8. 'Gaussian'
%                             9. 'Spiral'  
%       options        --- Options 
%              .para   --- parameter for some type. Detai look at each type
%
%
% Output: 
%       Data      --- preserving the older data set information and
%                      Diffusionmap parameter
%
%
%  Modification
%  WeiX Dec 5th 2014 create First Edition
%  WeiX Jul 16yh 2016 add Spiral
%
% Reference:
% 1) by Todd Wig1man, Department of Mathematics, University of Minnesota
%    E-mail wig1man@math.ucla.edu with comments & questions.
%    MANI Website: hg1p://www.math.ucla.edu/~wig1man/mani/index.html
% 
%
%% Initialization
switch nargin
    case 1
        error('Type of data not assigned');
    case 2        
        options.exit=0;
        warning('Default options is use');
    case 3 
        
    otherwise 
        error('Too much input')
end

       
     
            
switch type
    case 'SwissRoll'
        
        % Local parameters & Preprocess
        theta_min=3/2*pi;
        theta_max=9/2*pi;
        
        width_min=-10;
        width_max=10;
        
        length_arc=pi/2*(theta_max^2-theta_min^2); % This equation to calculte the perimeter of a spiral curve is weong! Use a simple approximation below.     
        length_arc=60;
        length_width=width_max-width_min;
        
%         Num_arc=round(length_arc/(length_arc*length_width)*N);
        Num_arc=round(sqrt(N/(length_arc*length_width))*length_arc);    
        Num_width=round(sqrt(N/(length_arc*length_width))*length_width);   
%         Num_width=round(length_width/(length_arc*length_width)*N);
        

        Num=Num_arc*Num_width;
        if (Num~=N)
            warning('Number of data point automatically adapts to %d',Num);
        end
        
        g1v=theta_min:(theta_max-theta_min)/(Num_arc-1):theta_max;   %grid variation dimension 1
        g2v=width_min:(width_max-width_min)/(Num_width-1):width_max; %grid variation dimension 2
        [g1,g2] = meshgrid(g1v, g2v);
        g1=g1(:);       % grid dimension 1
        g2=g2(:);       % grid dimension 2
        
%         % New version to generate g1,g2.! But Fail!
%         g1 = linspace(theta_min,theta_max,Num_arc);
%         g2 = linspace(width_min,width_max,Num_width);
        
        Data.X = [g1,g2];
        Data.Y = [g1.*cos(g1), g2,g1.*sin(g1)];
        Data.ColorVector = g1;
        Data.SizeVector = (g2-width_min+1)*5;

    case 'SwissHole'
        % Swiss Roll w/ hole example taken from Donoho & Grimes
        
        % Local parameters & Preprocess
        theta_min=3/2*pi;
        theta_max=9/2*pi;
        
        width_min=-10;
        width_max=10;
        
        length_arc=pi/2*(theta_max^2-theta_min^2);    
        length_arc=60;
        length_width=width_max-width_min;
        
%         Num_arc=round(length_arc/(length_arc*length_width)*N);
        Num_arc=round(sqrt(N/(length_arc*length_width))*length_arc);    
        Num_width=round(sqrt(N/(length_arc*length_width))*length_width);   
%         Num_width=round(length_width/(length_arc*length_width)*N);

        
        Num=Num_arc*Num_width;
        if (Num~=N)
            warning('Number of data point automatically adapts to %d',Num);
        end
        
        g1v=theta_min:(theta_max-theta_min)/(Num_arc-1):theta_max;   %grid variation dimension 1
        g2v=width_min:(width_max-width_min)/(Num_width-1):width_max; %grid variation dimension 2
        [g1,g2] = meshgrid(g1v, g2v);
        g1=g1(:);       % grid dimension 1
        g2=g2(:);       % grid dimension 2       
        
%         for i=1:Num_arc
%             for j=1:Num_width
%                 if (i>Num_arc/2)&(j<Num_width*3/4)&(j>Num_width*1/4)
%                     g1((i-1)*Num_arc+(j-1)*Num_width)=[];
%                     g2((i-1)*Num_arc+(j-1)*Num_width)=[];
%                 end
%             end
%         end
        index=[];
        for i=1:Num
           if (g1(i)>(theta_min+(theta_max-theta_min)/2)) & (g2(i)>(width_min+(width_max-width_min)*1/4)) & (g2(i)<(width_min+(width_max-width_min)*3/4))
               index=[index;i];           
           end
        end
        g1(index)=[];
        g2(index)=[];
        
        [Num,~]=size(g1);
        if (Num~=N)
            warning('Number of data point automatically adapts to %d',Num);
        end
        
        Data.X = [g1,g2];
        Data.Y = [g1.*cos(g1), g2, g1.*sin(g1)];     
        Data.ColorVector = g1+1e-5;   % Avoid zero value
%         Data.SizeVector = (g2-width_min+1)*5;
        Data.SizeVector = 20*ones(Num,1);
        
    case 'CornerPlanes'
        k = 1;
        xMax = floor(sqrt(N));
        yMax = ceil(N/xMax);
        cornerPoint = floor(yMax/2);
        for x = 0:xMax
            for y = 0:yMax
                if y <= cornerPoint
                    X(k,:) = [x,y,0];
                    ColorVector(k) = y;
                else
                    X(k,:) = [x,cornerPoint+(y-cornerPoint)*cos(pi*1/180),(y-cornerPoint)*sin(pi*1/180)];
                    ColorVector(k) = y;
                end;
                k = k+1;
            end;
        end;
        Data.Y = X;
        Data.ColorVector = ColorVector';
        Data.SizeVector = 10*ones((xMax+1)*(yMax+1),1);

    case 'PuncturedSphere'  %by Saul & Roweis
        inc = 9/sqrt(N);   %inc = 1/4;
        [xx,yy] = meshgrid(-5:inc:5);
        rr2 = xx(:).^2 + yy(:).^2;
        [tmp ii] = sort(rr2);
        Y = [xx(ii(1:N))'; yy(ii(1:N))'];
        a = 4./(4+sum(Y.^2));
        Data.X = Y';
        Data.Y = [a.*Y(1,:); a.*Y(2,:); 1*2*(1-a)]';
        Data.ColorVector = Data.Y(:,3);
        Data.SizeVector = 10*ones(N,1);
 
    case 'TwinPeaks'     % Twin Peaks by Saul & Roweis
        
        x_min=-1;
        x_max=1;
        
        y_min=-1;
        y_max=1;
        
        length_x=x_max-x_min;    
        length_y=y_max-y_min;  
        
        Num_x=round(sqrt(N/(length_x*length_y))*length_x);    
        Num_y=round(sqrt(N/(length_x*length_y))*length_y);   
        Num=Num_x*Num_y;  
        if (Num~=N)
            warning('Number of data point automatically adapts to %d',Num);
        end
        
        g1v=x_min:(x_max-x_min)/(Num_x-1):x_max;   %grid variation dimension 1
        g2v=y_min:(y_max-y_min)/(Num_y-1):y_max; %grid variation dimension 2
        [g1,g2] = meshgrid(g1v, g2v);
        g1=g1(:);       % grid dimension 1
        g2=g2(:);       % grid dimension 2
        
        Data.X = [g1,g2];
        Data.Y = [g1,g2,sin(pi*g1).*tanh(3*g2)];
%         Data.Y(:,3) = 1 * Data.Y(:,3);
        Data.ColorVector = Data.Y(:,3);
        Data.SizeVector = 10*ones(Num,1);
        
    case '3DClusters'       
%         options.ss = 8;
        if ~isfield(options,'para'), options.para = 3; end                 % Default Cluster number;
        numClusters = options.para;
        if (numClusters~=round(numClusters)), error('Cluster number must be integer');end
%         numClusters = max(1,numClusters);
        Centers = 10*rand(numClusters,3);
        D = L2_distance(Centers',Centers',1);
        minDistance = min(D(find(D>0)));
        k = 1;
        N2 = N - (numClusters-1)*9;
        for i = 1:numClusters
            for j = 1:ceil(N2/numClusters)
               X(k,1:3) = Centers(i,1:3)+(rand(1,3)-0.5)*minDistance/sqrt(12);
               ColorVector(k) = i;
               k = k + 1;
           end;
           % Connect clusters with straight line.
           if i < numClusters
               for t = 0.1:0.1:0.9
                    X(k,1:3) = Centers(i,1:3) + (Centers(i+1,1:3)-Centers(i,1:3))*t;
                    ColorVector(k) = 0;
                    k = k+1;
                end;
           end;
        end;
        Data.Y = X;
        Data.ColorVector = ColorVector;
        [Npoint,~]=size(Data.Y);
        Data.SizeVector = 10*ones(Npoint,1);

    case 'ToroidalHelix' % Toroidal Helix by Coifman & Lafon
        if ~isfield(options,'para'), options.para = 0.00; end                 % Default Cluster number;
        noiseSigma=options.para;   %noise parameter 
        if noiseSigma>0.1, warning('noise might be too large (>0.1)'); end
        t = (1:N)'/N;
        t = t.^(1)*2*pi;
        Data.X=t;
        Data.Y = [(2+cos(8*t)).*cos(t),(2+cos(8*t)).*sin(t),sin(8*t)]+noiseSigma*randn(N,3);
        Data.ColorVector = t;
        Data.SizeVector = 10*ones(N,1);

    case 'Gaussian'  % Gaussian randomly sampled
%         X = 1 * randn(N,3);
%         X(:,3) = 1 / (1^2 * 2 * pi) * exp ( (-X(:,1).^2 - X(:,2).^2) / (2*1^2) );
%         Data.Y = X;
%         Data.ColorVector = X(:,3);
%         Data.SizeVector = 10*ones(N,1);
        
        
        x_min=-1;
        x_max=1;
        
        y_min=-1;
        y_max=1;
        
        length_x=x_max-x_min;    
        length_y=y_max-y_min;  
        
        Num_x=round(sqrt(N/(length_x*length_y))*length_x);    
        Num_y=round(sqrt(N/(length_x*length_y))*length_y);   
        Num=Num_x*Num_y;  
        if (Num~=N)
            warning('Number of data point automatically adapts to %d',Num);
        end
        
        g1v=x_min:(x_max-x_min)/(Num_x-1):x_max;   %grid variation dimension 1
        g2v=y_min:(y_max-y_min)/(Num_y-1):y_max; %grid variation dimension 2
        [g1,g2] = meshgrid(g1v, g2v);
        g1=g1(:);       % grid dimension 1
        g2=g2(:);       % grid dimension 2
        
        Data.X=[g1,g2];
        Data.Y = [g1, g2, exp((-g1.^2 - g2.^2)/0.5) ];
%         Data.Y(:,3) = 1 * Data.Y(:,3);
        Data.ColorVector = Data.Y(:,3);
        Data.SizeVector = 10*ones(Num,1);
        
        
    case 'Spiral' 
        tt = (1:N)'/N;
        tt = tt.^(1)*6*pi;       
%         tt=rand(N,1)*6*pi;
        Data.X = [tt];
        Data.Y = [tt.*cos(tt),1*tt.*sin(tt),tt*10];
        Data.ColorVector = tt;
        Data.SizeVector = 10*ones(N,1);
        
        
        
    otherwise 
        error('No such type')
end

%     case 'OccludedDisks'  % Occluded disks
%         m = 20;   % Image size m x m.
%         Rd = 3;   % Disk radius.
%         Center = (m+1)/2;
%         u0 = zeros(m,m);
%         for i = 1:m
%             for j = 1:m
%                 if (Center - i)^2 + (Center - j)^2 < 1^2
%                     u0(i,j) = 1;
%                 end;
%             end;
%         end;
%         for diskNum = 1:N
%             DiskX(diskNum) = (m-1)*rand+1;
%             DiskY(diskNum) = (m-1)*rand+1;
%             u = u0;
%             for i = 1:m
%                 for j = 1:m
%                     if (DiskY(diskNum) - i)^2 + (DiskX(diskNum) - j)^2 < Rd^2
%                         u(i,j) = 1;
%                     end;
%                 end;
%             end;
%             X(diskNum,:) = reshape(u,1,m*m);
%         end;
%         % Since this is a special manifold, plot separately.
% %         axes(Data.maniAXES);
%         t = 0:0.1:2*pi+0.1;
% %         plot(1*cos(t)+Center,1*sin(t)+Center);
% %         axis([0.5 m+0.5 0.5 m+0.5]);
% %         hold on;
%         Data.ColorVector = (sqrt((DiskX-Center).^2+(DiskY-Center).^2))';
% %         scag1er(DiskX,DiskY,12,Data.ColorVector');
% %         hold off;
%         Data.Y = X;
%         [Npoint,~]=size(Data.Y);
%         Data.SizeVector = 10*ones(Npoint,1);
           

   
return



% --- L2_distance function
% Wrig1en by Roland Bunschoten, University of Amsterdam, 1999
function d = L2_distance(a,b,df)
if (size(a,1) == 1)
  a = [a; zeros(1,size(a,2))]; 
  b = [b; zeros(1,size(b,2))]; 
end
aa=sum(a.*a); bb=sum(b.*b); ab=a'*b; 
d = sqrt(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab);
d = real(d); 
if (df==1)
  d = d.*(1-eye(size(d)));
end
return 



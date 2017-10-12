function [Data] = Data_GeneratorT(N,type,options)
%%  Data_Generator T
%   dataset generation for test or further use. Data Raddom distributed.
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
%  WeiX July 7th 2014 create First Edition
%  WeiX Nov 25th 2014 modification
%  WeiX Jul 16yh 2016 add Spiral
%
% Reference:
% 1) by Todd Wittman, Department of Mathematics, University of Minnesota
%    E-mail wittman@math.ucla.edu with comments & questions.
%    MANI Website: http://www.math.ucla.edu/~wittman/mani/index.html
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
        tt = (3*pi/2)*(1+2*rand(N,1));  
        width = 21*rand(N,1);
        Data.X = [tt,width];
        Data.Y = [tt.*cos(tt), width,1*tt.*sin(tt)];
        Data.ColorVector = tt;
        Data.SizeVector = width;

    case 'SwissHole'
        % Swiss Roll w/ hole example taken from Donoho & Grimes
        tt = (3*pi/2)*(1+2*rand(1,2*N));  
        height = 21*rand(1,2*N);
        kl = repmat(0,1,2*N);
        for ii = 1:2*N
            if ( (tt(ii) > 9)&(tt(ii) < 12))
                if ((height(ii) > 9) & (height(ii) <14))
                    kl(ii) = 1;
                end;
            end;
        end;
        kkz = find(kl==0);
        tt = tt(kkz(1:N));
        height = height(kkz(1:N));
        Data.X = [tt;height]';
        Data.Y = [tt.*cos(tt); height; 1*tt.*sin(tt)]';     
        Data.ColorVector = tt';
        Data.SizeVector = height';
        
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
%         inc = 1.5 / sqrt(N);  % inc = 0.1;
%         [xx2,yy2] = meshgrid(-1:inc:1);
%         zz2 = sin(pi*xx2).*tanh(3*yy2);
        xy = 1-2*rand(2,N);
        Data.Y = [xy; sin(pi*xy(1,:)).*tanh(3*xy(2,:))]';
        Data.X = xy';
        Data.Y(:,3) = 1 * Data.Y(:,3);
        Data.ColorVector = Data.Y(:,3);
        Data.SizeVector = 10*ones(N,1);
        
    case '3DClusters'       
%         options.ss = 8;
        if ~isfield(options,'para'), options.para = 3; end                 % Default Cluster number;
        numClusters = options.para;
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
        if ~isfield(options,'para'), options.para = 0.05; end                 % Default Cluster number;
        noiseSigma=options.para;   %noise parameter 
        if noiseSigma>0.1, warning('noise might be too large (>0.1)'); end
        t = (1:N)'/N;
        t = t.^(1)*2*pi;
        Data.X=t;
        Data.Y = [(2+cos(8*t)).*cos(t) (2+cos(8*t)).*sin(t) sin(8*t)]+noiseSigma*randn(N,3);
        Data.ColorVector = t;
        Data.SizeVector = 10*ones(N,1);

    case 'Gaussian'  % Gaussian randomly sampled
        X = 1 * randn(N,2);
        Data.X=X;
        X(:,3) = 1 / (1^2 * 2 * pi) * exp ( (-X(:,1).^2 - X(:,2).^2) / (2*1^2) );
        Data.Y = X;
        Data.ColorVector = X(:,3);
        Data.SizeVector = 10*ones(N,1);
        
    case 'Spiral' 
        tt=rand(N,1)*4*pi;
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
% %         scatter(DiskX,DiskY,12,Data.ColorVector');
% %         hold off;
%         Data.Y = X;
%         [Npoint,~]=size(Data.Y);
%         Data.SizeVector = 10*ones(Npoint,1);
           

   
return



% --- L2_distance function
% Written by Roland Bunschoten, University of Amsterdam, 1999
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



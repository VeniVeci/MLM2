function [ ] = DisplayUV( U,V,P,Re,tf,dt,lx,ly,nx,ny )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

nt = ceil(tf/dt); dt = tf/nt;
x = linspace(0,lx,nx+1); hx = lx/nx;
y = linspace(0,ly,ny+1); hy = ly/ny;
[X,Y] = meshgrid(y,x);
%-----------------------------------------------------------------------
% initial conditions
% U = zeros(nx-1,ny); V = zeros(nx,ny-1);
% U_old = zeros(nx-1,ny); V_old = zeros(nx,ny-1);

% boundary conditions
uN = x*0;    vN = avg(x)*0;
uS = x*0;      vS = avg(x)*0;
uW = avg(y)*0; vW = y*0;
uE = avg(y)*0; vE = y*0;
%-----------------------------------------------------------------------
Ubc = dt/Re*([2*uS(2:end-1)' zeros(nx-1,ny-2) 2*uN(2:end-1)']/hx^2+...
      [uW;zeros(nx-3,ny);uE]/hy^2);
Vbc = dt/Re*([vS' zeros(nx,ny-3) vN']/hx^2+...
      [2*vW(2:end-1);zeros(nx-2,ny-1);2*vE(2:end-1)]/hy^2);

% fprintf('initialization')
Lp = kron(speye(ny),K1(nx,hx,1))+kron(K1(ny,hy,1),speye(nx));
Lp(1,1) = 3/2*Lp(1,1);
perp = symamd(Lp); Rp = chol(Lp(perp,perp)); Rpt = Rp';
Lu = speye((nx-1)*ny)+dt/Re*(kron(speye(ny),K1(nx-1,hx,2))+...
     kron(K1(ny,hy,3),speye(nx-1)));
peru = symamd(Lu); Ru = chol(Lu(peru,peru)); Rut = Ru';
Lv = speye(nx*(ny-1))+dt/Re*(kron(speye(ny-1),K1(nx,hx,3))+...
     kron(K1(ny-1,hy,2),speye(nx)));
perv = symamd(Lv); Rv = chol(Lv(perv,perv)); Rvt = Rv';
Lq = kron(speye(ny-1),K1(nx-1,hx,2))+kron(K1(ny-1,hy,2),speye(nx-1));
perq = symamd(Lq); Rq = chol(Lq(perq,perq)); Rqt = Rq';



rhs = reshape(diff([uW;U;uE])/hx+diff([vS' V vN']')'/hy,[],1);
p(perp) = -Rp\(Rpt\rhs(perp));
% P = reshape(p,nx,ny);


% stream function
rhs = reshape(diff(U')'/hy-diff(V)/hx,[],1);
q(perq) = Rq\(Rqt\rhs(perq));
Q = zeros(nx+1,ny+1);
Q(2:end-1,2:end-1) = reshape(q,nx-1,ny-1);
clf, contourf(avg(x),avg(y),P',20,'w-'), hold on
contour(x,y,Q',20,'k-');
Ue = [uS' avg([uW;U;uE]')' uN'];
Ve = [vW;avg([vS' V vN']);vE];
Len = sqrt(Ue.^2+Ve.^2+eps);
quiver(x,y,(Ue./Len)',(Ve./Len)',.4,'k-')
hold off, axis equal, axis([0 lx 0 ly])
p = sort(p); 

% caxis(p([8 end-7]))
caxis auto

title(sprintf('Re = %0.1g'))
colorbar
drawnow










return

function B = avg(A,k)
if nargin<2, k = 1; end
if size(A,1)==1, A = A'; end
if k<2, B = (A(2:end,:)+A(1:end-1,:))/2; else, B = avg(A,k-1); end
if size(A,2)==1, B = B'; end
return

function A = K1(n,h,a11)
% a11: Neumann=1, Dirichlet=2, Dirichlet mid=3;
A = spdiags([-1 a11 0;ones(n-2,1)*[-1 2 -1];0 a11 -1],-1:1,n,n)'/h^2;
return

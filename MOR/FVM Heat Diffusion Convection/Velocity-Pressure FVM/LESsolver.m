function x=LESsolver(A,b,option)
% Linear Equation System solver by specified some methods.
%    
% Options are:
% NGE   for ---Navie Gaussian Elimination solver----------
% TDMA  for ---TDMA solver--------------------------------  
% BSOR  for ---bandSOR solver-----------------------------
% GE    for ---Gaussian elimination with back substitution
% JCB   for ---Jacobi Method------------------------------
% GS    for ---Gauss-Seidel iterations--------------------
% SOR   for ---SOR iterations-----------------------------
% CG    for ---conjugate gradient method------------------
% NT    for ---Newton's method ---------------------------
% BD    for ---Broyden's method ---------------------------
% 
% 
% Modifications:
% 90-May-2014, WeiX, first edition 

[rows,~]=size(b);

switch option
    
   case 'NGE'
       x=naiv_gauss(A,b);
   case 'TDMA'
       x=TDMAsolver(A,b);
   case 'BSOR'
       x0=rand(rows,1);
       omega=.5;
       epsilon=0.001;
       [x,Iteration,Error] = bandSOR(A,x0,b,omega,epsilon);   
   case 'GE'
        x = gauss_elim ( A, b );
   case 'JCB'
        xold=zeros(rows,1); 
        TOL=1e-6;
        Nmax=1e4;
        x = jacobi ( A, b, xold, TOL, Nmax );
    case 'GS'
        xold=zeros(rows,1); 
        TOL=1e-6;
        Nmax=1e4;
      x = gauss_seidel ( A, b, xold, TOL, Nmax );
    case 'SOR'
        xold=zeros(rows,1); 
        omega=1;
        TOL=1e-6;
        Nmax=1e4;
        x = sor ( A, b, xold, omega, TOL, Nmax );
    case 'CG'
        x=zeros(rows,1); 
        TOL=1e-6;
        Nmax=1e4;        
        x = conj_grad ( A, b, x, TOL, Nmax );
    case 'NT'
    case 'BD'
   otherwise
      disp('No such method');
end



end




%% ----------Navie Gaussian Elimination solver------------------------
function x = naiv_gauss(A,b);
n = length(b); x = zeros(n,1);
for k=1:n-1 % forward elimination
      for i=k+1:n
           xmult = A(i,k)/A(k,k);
           for j=k+1:n
             A(i,j) = A(i,j)-xmult*A(k,j);
           end
           b(i) = b(i)-xmult*b(k);
      end
end
 % back substitution
x(n) = b(n)/A(n,n);
for i=n-1:-1:1
   sum = b(i);
   for j=i+1:n
     sum = sum-A(i,j)*x(j);
   end
   x(i) = sum/A(i,i);
end
end

%% -----------------------------TDMA solver--------------------------------
function X=TDMAsolver(A,b)

m=length(b);                 % m is the number of rows
X=zeros(m,1);
A(1,2)= A(1,2)  ./ A(1,1);    % Division by zero risk.
b(1)=  b(1)    ./ A(1,1);    % Division by zero would imply a singular matrix

for i=2:m-1
    temp=  A(i,i) - A(i,i-1) .* A(i-1,i);
    A(i,i+1)=  A(i,i+1)  ./ temp;
    b(i)= ( b(i) - A(i,i-1) .* b(i-1) )  ./ temp;
end 

i=m;
X(m)=(b(i) - A(i,i-1) .* b(i-1))  ./ (A(i,i) - A(i,i-1) .* A(i-1,i));

for i=m-1:-1:1
X(i)=  -A(i,i+1) .* X(i+1) + b(i);
end
end




%% --------------------------bandSOR solver--------------------------------
function [x,Iteration,Error] = bandSOR(A,x0,b,omega,epsilon)
% Goal: This code solves the linear system Ax=b, where A is a symmetric banded
%      matrix, using banded SOR. A is first stored in compact storage mode
%      and next SOR is applied to the compact storage.
% Input: This code accepts a symmetric banded matrix A, initial guess x0, 
%       vector b, SOR parameter omega and tolerance epsilon.
% Output: This code outputs the SOR-computed solution x, number of iteration
%        called Iteration and error called Error, which is the Euclidean
%        distance between kth iterate and k+1 iterate just when Error <
%        epsilon.
% Author: Michael Akinwumi
%  (c) 2011
Error=1;
Iteration=0;
x=x0;
x=x(:);b=b(:);
B=compactstorage(A);
dim=size(B);
while Error>epsilon
    Iteration=Iteration+1;
    temp = x;
     for i=1:length(x)
        before = 0;
        for j=1:(i-1)
            if dim(2)>=i+1-j 
                before = before + B(j,i+1-j)*x(j);
            end;
        end;
        after = 0;
        for j=(i+1):length(x)
            if dim(2)>=j+1-i 
                after = after + B(i,j+1-i)*x(j);
            end;
        end
        x(i) = (omega/B(i,1))*(b(i)- before - after) + (1-omega)*x(i);
    end;
Error = norm(x-temp);
end;
end

function B = compactstorage(A)
% This function stores symmetric banded matrix A in a 
% compact form bA in such a way that only the main diagonal,
% and the nonzero superdiagonals are stored. The first column 
% of bA corresponds to the main diagonal of A and the subsequent 
% columns of bA correspond to superdiagonals of A.
% Input:upper or lower bandwidth p and a symmetric matrix A
% Output: bA, compressed form of A
dim=size(A);
if ~(dim(1)==dim(2))
    error('A must be square')
end
if (all((all(A)~=all(A'))))
    error('A must be symmetric')
end
if ~(all(eig(A))> 0)
    error('Matrix is at least not positive definite')
end
c=find(A(1,1:dim(1))~=0);
B=zeros(dim(1),c(end));
n=dim(1);p=c(end)-1;
for i=1:n
if i<=n-p
for j=i:p+i
B(i,j-i+1)=A(i,j);
end
else 
for j=i:n
B(i,j-i+1)=A(i,j);
end
end
end
end

%% ------------Gaussian elimination with back substitution----------------
function x = gauss_elim ( A, b )
%GAUSS_ELIM   solve the linear system Ax = b using Gaussian elimination
%             with back substitution
%
%     calling sequences:
%             x = gauss_elim ( A, b )
%             gauss_elim ( A, b )
%
%     inputs:
%             A       coefficient matrix for linear system
%                     (matrix must be square)
%             b       right-hand side vector
%
%     output:
%             x       solution vector (i.e., vector for which Ax = b)
%
%     NOTE:
%             this is intended as a demonstration routine - no pivoting
%             strategy is implemented to reduce the effects of roundoff
%             error
%

    [nrow ncol] = size ( A );
    if ( nrow ~= ncol )
       disp ( 'gauss_elim error: Square coefficient matrix required' );
       return;
    end;
    nb = length ( b );
    if ( nrow ~= nb )
       disp ( 'gauss_elim error: Size of b-vector not compatible with matrix dimension' )
       return;
    end;

x = zeros ( 1, nrow );

%
%    Gaussian elimination
%

    for i = 1 : nrow - 1
        if ( A(i,i) == 0 )
           t =  min ( find ( A(i+1:nrow,i) ~= 0 ) + i );
           if ( isempty(t) ) 
              disp ( 'gauss_elim error: A matrix is singular' );
              return
           end;
           temp = A(i,:);     tb = b(i);
           A(i,:) = A(t,:);   b(i) = b(t);
           A(t,:) = temp;     b(t) = tb;
        end;
        for j = i+1 : nrow
            m = -A(j,i) / A(i,i);
            A(j,i) = 0;
            A(j, i+1:nrow) = A(j, i+1:nrow) + m * A(i, i+1:nrow);
            b(j) = b(j) + m * b(i);
        end;
    end;

%
%    back substitution
%

    x(nrow) = b(nrow) / A(nrow, nrow);
    for i = nrow - 1 : -1 : 1
        x(i) = ( b(i) - sum ( x(i+1:nrow) .* A(i, i+1:nrow) ) ) / A(i,i);
    end;

end

%% -------------------Jacobi Method---------------------------------------
function x = jacobi ( A, b, xold, TOL, Nmax )

%JACOBI       approximate the solution of the linear system Ax = b by applying
%             the Jacobi method (simultaneous relaxation)
%
%     calling sequences:
%             x = jacobi ( A, b, xold, TOL, Nmax )
%             jacobi ( A, b, xold, TOL, Nmax )
%
%     inputs:
%             A       coefficient matrix for linear system - must be a
%                     square matrix
%             b       right-hand side vector for linear system
%             xold    vector containing initial guess for solution of 
%                     linear system
%             TOL     convergence tolerance - applied to maximum norm of
%                     difference between successive approximations
%             NMax    maximum number of iterations to be performed
%
%     output:
%             x       approximate solution of linear system
%
%     NOTE:
%             if JACOBI is called with no output arguments, the 
%             iteration number and the current approximation to the 
%             solution are displayed
%
%             if the maximum number of iterations is exceeded, a meesage
%             to this effect will be displayed and the current approximation 
%             will be returned in the output value
%

n = length ( b );
[r c] = size ( A );
if ( c ~= n ) 
   disp ( 'jacobi error: matrix dimensions and vector dimension not compatible' )
   return
end;
xnew = zeros ( 1, n );

if ( nargout == 0 )
   s = sprintf ( '%3d \t %10f ', 0, xold(1) );
   for j = 2 : n 
	   s = sprintf ( '%s%10f ', s, xold(j) );
   end;
   disp ( s );
end;

for its = 1 : Nmax
    xnew(1) = ( b(1) - sum ( A(1,2:n) .* xold(2:n) ) ) / A(1,1);
	for i = 2 : n-1
	    xnew(i) = ( b(i) - sum ( A(i,1:i-1) .* xold(1:i-1) ) - sum ( A(i,i+1:n) .* xold(i+1:n) ) ) / A(i,i);
	end;
	xnew(n) = ( b(n) - sum ( A(n,1:n-1) .* xold(1:n-1) ) ) / A(n,n);
	
    if ( nargout == 0 )
	   s = sprintf ( '%3d \t %10f ', its, xnew(1) );
	   for j = 2 : n 
	       s = sprintf ( '%s%10f ', s, xnew(j) );
	   end;
	   disp ( s );
	end;

    conv = max ( abs ( xnew - xold ) );
	if ( conv < TOL )
	   x = xnew;
	   return
	else
	   xold = xnew;
	end;
end;
disp ( 'jacobi error: maximum number of iterations exceeded' );
if ( nargout == 1 ) x = xnew; end;
	   
end

%% --------------- Gauss-Seidel iterations---------------------------------
function x = gauss_seidel ( A, b, xold, TOL, Nmax )

%GAUSS_SEIDEL approximate the solution of the linear system Ax = b by applying
%             the Gauss-Seidel method (successive relaxation)
%
%     calling sequences:
%             x = gauss_seidel ( A, b, xold, TOL, Nmax )
%             gauss_seidel ( A, b, xold, TOL, Nmax )
%
%     inputs:
%             A       coefficient matrix for linear system - must be a
%                     square matrix
%             b       right-hand side vector for linear system
%             xold    vector containing initial guess for solution of 
%                     linear system
%             TOL     convergence tolerance - applied to maximum norm of
%                     difference between successive approximations
%             NMax    maximum number of iterations to be performed
%
%     output:
%             x       approximate solution of linear system
%
%     NOTE:
%             if GAUSS_SEIDEL is called with no output arguments, the 
%             iteration number and the current approximation to the 
%             solution are displayed
%
%             if the maximum number of iterations is exceeded, a meesage
%             to this effect will be displayed and the current approximation 
%             will be returned in the output value
%

n = length ( b );
[r c] = size ( A );
if ( c ~= n ) 
   disp ( 'gauss_seidel error: matrix dimensions and vector dimension not compatible' )
   return
end;
xnew = zeros ( 1, n );

if ( nargout == 0 )
   s = sprintf ( '%3d \t %10f ', 0, xold(1) );
   for j = 2 : n 
	   s = sprintf ( '%s%10f ', s, xold(j) );
   end;
   disp ( s );
end;

for its = 1 : Nmax
    xnew(1) = ( b(1) - sum ( A(1,2:n) .* xold(2:n) ) ) / A(1,1);
	for i = 2 : n-1
	    xnew(i) = ( b(i) - sum ( A(i,1:i-1) .* xnew(1:i-1) ) - sum ( A(i,i+1:n) .* xold(i+1:n) ) ) / A(i,i);
	end;
	xnew(n) = ( b(n) - sum ( A(n,1:n-1) .* xnew(1:n-1) ) ) / A(n,n);
	
    if ( nargout == 0 )
	   s = sprintf ( '%3d \t %10f ', its, xnew(1) );
	   for j = 2 : n 
	       s = sprintf ( '%s%10f ', s, xnew(j) );
	   end;
	   disp ( s );
	end;

    conv = max ( abs ( xnew - xold ) );
	if ( conv < TOL )
	   x = xnew;
	   return
	else
	   xold = xnew;
	end;
end;
disp ( 'gauss_seidel error: maximum number of iterations exceeded' );
if ( nargout == 1 ) x = xnew; end;
	   
end

%% ---------------SOR iterations------------------------------------------
function x = sor ( A, b, xold, omega, TOL, Nmax )

%SOR          approximate the solution of the linear system Ax = b by applying
%             the SOR method (successive over-relaxation)
%
%     calling sequences:
%             x = sor ( A, b, xold, omega, TOL, Nmax )
%             sor ( A, b, xold, omega, TOL, Nmax )
%
%     inputs:
%             A       coefficient matrix for linear system - must be a
%                     square matrix
%             b       right-hand side vector for linear system
%             xold    vector containing initial guess for solution of 
%                     linear system
%             omega   relaxation parameter
%             TOL     convergence tolerance - applied to maximum norm of
%                     difference between successive approximations
%             NMax    maximum number of iterations to be performed
%
%     output:
%             x       approximate solution of linear system
%
%     NOTE:
%             if SOR is called with no output arguments, the 
%             iteration number and the current approximation to the 
%             solution are displayed
%
%             if the maximum number of iterations is exceeded, a meesage
%             to this effect will be displayed and the current approximation 
%             will be returned in the output value
%

n = length ( b );
[r c] = size ( A );
if ( c ~= n ) 
   disp ( 'sor error: matrix dimensions and vector dimension not compatible' )
   return
end;
xnew = zeros ( 1, n );

if ( nargout == 0 )
   s = sprintf ( '%3d \t %10f ', 0, xold(1) );
   for j = 2 : n 
	   s = sprintf ( '%s%10f ', s, xold(j) );
   end;
   disp ( s );
end;

for its = 1 : Nmax
    xnew(1) = ( 1 - omega ) * xold(1) + omega * ( b(1) - sum ( A(1,2:n) .* xold(2:n) ) ) / A(1,1);
	for i = 2 : n-1
	    xnew(i) = ( 1 - omega ) * xold(i) + omega * ( b(i) - sum ( A(i,1:i-1) .* xnew(1:i-1) ) - sum ( A(i,i+1:n) .* xold(i+1:n) ) ) / A(i,i);
	end;
	xnew(n) = ( 1 - omega ) * xold(n) + omega * ( b(n) - sum ( A(n,1:n-1) .* xnew(1:n-1) ) ) / A(n,n);
	
    if ( nargout == 0 )
	   s = sprintf ( '%3d \t %10f ', its, xnew(1) );
	   for j = 2 : n 
	       s = sprintf ( '%s%10f ', s, xnew(j) );
	   end;
	   disp ( s );
	end;

    conv = max ( abs ( xnew - xold ) );
	if ( conv < TOL )
	   x = xnew;
	   return
	else
	   xold = xnew;
	end;
end;
disp ( 'sor error: maximum number of iterations exceeded' );
if ( nargout == 1 ) x = xnew; end;
end

%% ---------------conjugate gradient method--------------------------------
function x = conj_grad ( A, b, x, TOL, Nmax )

%CONJ_GRAD    approximate the solution of the linear system Ax = b by applying
%             the conjugate gradient method
%
%     calling sequences:
%             x = conj_grad ( A, b, x, TOL, Nmax )
%             conj_grad ( A, b, x, TOL, Nmax )
%
%     inputs:
%             A       coefficient matrix for linear system - must be a
%                     square matrix
%             b       right-hand side vector for linear system - must be
%                     a column vector
%             x       column vector containing initial guess for solution of 
%                     linear system
%             TOL     convergence tolerance - applied to Euclidean norm of
%                     residual vector
%             NMax    maximum number of iterations to be performed
%
%     output:
%             x       approximate solution of linear system
%
%     NOTE:
%             if CONJ_GRAD is called with no output arguments, the 
%             iteration number and the current approximation to the 
%             solution are displayed
%
%             if the maximum number of iterations is exceeded, a meesage
%             to this effect will be displayed and the current approximation 
%             will be returned in the output value
%

n = length ( b );
[r c] = size ( A );
if ( c ~= n ) 
   disp ( 'conj_grad error: matrix dimensions and vector dimension not compatible' )
   return
end;

if ( nargout == 0 )
   s = sprintf ( '%3d \t %10f ', 0, x(1) );
   for j = 2 : n 
	   s = sprintf ( '%s%10f ', s, x(j) );
   end;
   disp ( s );
end;
r = A * x - b;
delta0 = r' * r;
d = -r;

for its = 1 : Nmax
    h = A * d;
	lambda = delta0 / ( d' * h );
	x = x + lambda * d;
	r = r + lambda * h;
	delta1 = r' * r;
	
    if ( nargout == 0 )
	   s = sprintf ( '%3d \t %10f ', its, x(1) );
	   for j = 2 : n 
	       s = sprintf ( '%s%10f ', s, x(j) );
	   end;
	   disp ( s );
	end;

	if ( sqrt(delta1) < TOL )
	   return
	else
	   alpha = delta1 / delta0;
	   delta0 = delta1;
	   d = -r + alpha * d;
	end;
end;
disp ( 'conj_grad error: maximum number of iterations exceeded' );
end

%% ----- Newton's method --------------------------------------------
function y = newton_sys ( F, J, x0, TOL, Nmax )

%NEWTON_SYS   solve the system of nonlinear equations F(x) = 0 using 
%             Newton's method
%
%     calling sequences:
%             y = newton_sys ( F, J, x0, TOL, Nmax )
%             newton_sys ( F, J, x0, TOL, Nmax )
%
%     inputs:
%             F       vector-valued function of a vector argument which
%                     defines the system of equations to be solved
%             J       matrix-valued function which computes the Jacobian 
%                     associated with the function F
%             x0      vector containing initial guess for solution of 
%                     nonlinear system
%             TOL     convergence tolerance - applied to maximum norm of
%                     difference between successive approximations
%             NMax    maximum number of iterations to be performed
%
%     output:
%             y       approximate solution of nonlinear system
%
%     dependencies:
%             this routine uses both LUfactor and LUsolve
%
%     NOTE:
%             if NEWTON_SYS is called with no output arguments, each 
%             approximation to the solution is displayed
%
%             if the maximum number of iterations is exceeded, a meesage
%             to this effect will be displayed and the current approximation 
%             will be returned in the output value
%

old = x0;
for i = 1 : Nmax
    Fold = feval(F,old);
	Jold = feval(J,old);
	[lu pvt] = LUfactor ( Jold );
	dx = LUsolve ( lu, -Fold, pvt );
    new = old + dx;
	
	if ( nargout == 0 )
	   disp ( new )
	end
	
	if ( max(abs(dx)) < TOL ) 
	   if ( nargout == 1 )
	      y = new;
	   end
	   return
	else
	   old = new;
	end
end

disp('newton_sys error: Maximum number of iterations exceeded');
if ( nargout == 1 ) y = new; end;
end

%% -----------------Broyden's method --------------------------------------
function y = broyden ( F, J, x0, TOL, Nmax )

%BROYDEN      solve the system of nonlinear equations F(x) = 0 using 
%             Broyden's method
%
%     calling sequences:
%             y = broyden ( F, J, x0, TOL, Nmax )
%             broyden ( F, J, x0, TOL, Nmax )
%
%     inputs:
%             F       vector-valued function of a vector argument which
%                     defines the system of equations to be solved
%             J       matrix-valued function which computes the Jacobian 
%                     associated with the function F
%             x0      vector containing initial guess for solution of 
%                     nonlinear system
%             TOL     convergence tolerance - applied to maximum norm of
%                     difference between successive approximations
%             NMax    maximum number of iterations to be performed
%
%     output:
%             y       approximate solution of nonlinear system
%
%
%     NOTE:
%             if BROYDEN is called with no output arguments, each 
%             approximation to the solution is displayed
%
%             if the maximum number of iterations is exceeded, a meesage
%             to this effect will be displayed and the current approximation 
%             will be returned in the output value
%

Fold = feval(F,x0)';
Jold = feval(J,x0);
A0 = inv ( Jold );
dx = -A0 * Fold;
x0  = x0 + dx;
if ( nargout == 0 )
	disp ( x0' )
end

for i = 2 : Nmax
    Fnew = feval(F,x0)';
	dy = Fnew - Fold;
	u = A0 * dy;
	v = dx' * A0;
	denom = dx' * u;
	A0 = A0 + ( dx - u ) * v / denom;
	dx = -A0 * Fnew;
    x0 = x0 + dx;
	
	if ( nargout == 0 )
	   disp ( x0' )
	end
	
	if ( max(abs(dx)) < TOL ) 
	   if ( nargout == 1 )
	      y = x0;
	   end
	   return
	else
	   Fold = Fnew;
	end
	
end

disp('broyden error: Maximum number of iterations exceeded');
if ( nargout == 1 ) y = x0; end;
end




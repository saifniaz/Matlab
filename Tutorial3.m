function Tutorial3
format short

f = @(x) 2*x.^3 + 4*x.^2 - 4*x - 6;
exact = fzero(f,1)
%[br,vecbr] = BM(f,1,2,15)
%g = @(x)  (x^3+80)/(2*x^2);
g = @(x)(exp(x)-sin(x))/3;
%[bfp, vecbfp] = FP(g,0,5,1e-5)
[p1, p2, p, P] = MM(f,1.1,1.3,1.5,1e-5,1e-5,3)
%[bsm, vecbsm] = SM(f,0,6,1e-5,5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xc,x] = BM(f, a, b,Maxit)
if (sign(f(a))*sign(f(b)) > 0);
error('BM condition not satisfied.')
return
end

 Niter_BM = 0;
for k = 1 : Maxit
c = (a+b)/2;
x(k) =c;
%cond = abs(f(c)); % a number to compare with tolerance
if f(c) == 0; xc = c; break; end % solution is found
if sign(f(c))*sign(f(b))< 0; a = c; else b = c; end % get new interval contaning solution
end
xc = (a + b)/2; %new midpoint is best estimate
x(k+1)=xc;
 %Niter_BM = Niter_BM+1
end
%% 2. Fixed-Point Iteration Method
% Program 1.2 Fixed-Point Iteration
% Computes approximate solution of g(x) = x
% Input: function handle g, starting guess x0,
% number of iteration steps k
% Output: Approximate solution xc
function [xc,x] = FP(g, x0,Maxit, tol)
     Niter_FP = 0;
x(1) = x0;
for i = 1: Maxit
x(i+1) = g(x(i));
%Niter_FP = Niter_FP+1;
%if (N > ceil(k/2) && abs(x(i+1)-x(i))>1); error('not converging'); break; end
if abs(x(i+1)-x(i))<tol; break; end
end
xc = x(i+1);
% Niter_FP = Niter_FP+1
end
function [xs,x] = SM(f,a,b,tol, Maxit)
    x(1) = a; x(2) = b;
    err=100;
    Niter_SM = 0;
    for k = 2 : Maxit
        x(k+1) = x(k) - (f(x(k))*(x(k)-x(k-1)))/(f(x(k))-f(x(k-1)));
    
        if abs(x(k+1)-x(k))<tol; break; end
        
      %  Niter_SM = Niter_SM+1
    end
    xs = x(k+1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p2,y2,err,P] = MM(f,p0,p1,p2,delta,epsilon,max1)
%---------------------------------------------------------------------------
%MULLER   Muller's method is used to locate a root
% Sample calls
%   [p2,y2,err] = muller('f',p0,p1,p2,delta,epsilon,max1)
%   [p2,y2,err,P] = muller('f',p0,p1,p2,delta,epsilon,max1)
% Inputs
%   f         name of the function
%   p0        starting value
%   p1        starting value
%   p2        starting value
%   delta     convergence tolerance for p2
%   epsilon   convergence tolerance for y2
%   max1      maximum number of iterations
% Return
%   p2        solution: the root
%   y2        solution: the function value
%   err       error estimate in the solution for p2
%   P         History vector of the iterations
%
% NUMERICAL METHODS: MATLAB Programs, (c) John H. Mathews 1995
% To accompany the text:
% NUMERICAL METHODS for Mathematics, Science and Engineering, 2nd Ed, 1992
% Prentice Hall, Englewood Cliffs, New Jersey, 07632, U.S.A.
% Prentice Hall, Inc.; USA, Canada, Mexico ISBN 0-13-624990-6
% Prentice Hall, International Editions:   ISBN 0-13-625047-5
% This free software is compliments of the author.
% E-mail address:      in%"mathews@fullerton.edu"
%
% Algorithm 2.8 (Muller's Method).
% Section	2.5, Aitken's Process & Steffensen's & Muller's Methods, Page 97
%---------------------------------------------------------------------------

P(1) = p0;
P(2) = p1;
P(3) = p2;
y0  = feval(f,p0);
y1  = feval(f,p1);
y2  = feval(f,p2);
for k=1:max1,
  h0 = p0 - p2;
  h1 = p1 - p2;
  c  = y2;
  e0 = y0 - c;
  e1 = y1 - c;
  det1 = h0*h1*(h0-h1);
  a  = (e0*h1 - h0*e1)/det1;
  b  = (h0^2*e1 - h1^2*e0)/det1;
  if  b^2 > 4*a*c,
    disc = sqrt(b^2 - 4*a*c);
  else
    disc = 0;
  end
  if b < 0, disc = - disc; end
  z = - 2*c/(b + disc);
  p3 = p2 + z;
  if  abs(p3-p1) < abs(p3-p0),
    u = p1;
    p1 = p0;
    p0 = u;
    v = y1;
    y1 = y0;
    y0 = v;
  end
  if  abs(p3-p2) < abs(p3-p1),
    u = p2;
    p2 = p1;
    p1 = u;
    v = y2;
    y2 = y1;
    y1 = v;   
  end
  p2 = p3;
  y2 = feval(f,p2);
  P = [P,p2];
  err = abs(z);
  relerr = err/(abs(p3)+eps);
  if (err<delta)|(relerr<delta)|(abs(y1)<epsilon), break, end
end
end
end

function Tutorial2
%%
% function definitions using function handle
f = @(x)  x.^6 - x - 1;
g = @(x) x-x^3;
% Interval 
a = 1;
b = 2;
% Initial guess
x0 = -1;
% Tolerance
Tol = 1e-8;
% Function Calling
%BMroot = bisect(f,a,b,Tol)
Fzeroot = fzero(f,x0)
%FProot = fpi(g,x0,100)
%%
% The Wilkinsons Polynomials
%web('http://blogs.mathworks.com/cleve/2013/03/04/wilkinsons-polynomials/')

%% Program 1.1 Bisection Method
% Computes approximate solution of f(x)=0
% Input: function handle f; a,b such that f(a) ? f(b) < 0,
% and tolerance tol
% Output: Approximate solution xc
function xc = bisect(f, a, b,tol)
if (sign(f(a))*sign(f(b)) > 0);
error('BM condition not satisfied.');
return
end
cond = 100;
while cond > tol
c = (a+b)/2;
cond = abs(f(c)); % a number to compare with tolerance
if f(c) == 0; xc = c; break; end % solution is found
if sign(f(c))*sign(f(b))< 0; a = c; else b = c; end % get new interval contaning solution
end
xc = (a + b)/2; %new midpoint is best estimate
end
%% 2. Fixed-Point Iteration Method
% Program 1.2 Fixed-Point Iteration
% Computes approximate solution of g(x) = x
% Input: function handle g, starting guess x0,
% number of iteration steps k
% Output: Approximate solution xc
function xc1 = fpi(g, x0, k)
    N=1;
x(1) = x0;
for i = 1: k
x(i+1) = g(x(i));
N=N+1;
if (N > ceil(k/2) && abs(x(i+1)-x(i))>1); error('not converging'); break; end
end
xc1 = x(i+1);
end
end

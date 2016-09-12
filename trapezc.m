function It = trapezc(a,b,m,f)
% It = trapezc(a,b,m,f)
% INPUT:
% a  == left endpoint of interval
% b  == right endpoint of interval
% m  == number of subintervals
% f  == function handle for integrand (vectorised)
% OUTPUT:
% It == computed approximation of integral
h = (b-a)/m;
x = [a:h:b];
dim = max(size(x));
y = f(x); 
if size(y)==1, y=diag(ones(dim))*y; end; 
It = h*(0.5*y(1)+sum(y(2:m))+0.5*y(m+1));

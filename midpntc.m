function Im = midpntc(a,b,m,f)
% Im = midpntc(a,b,m,f)
% INPUT:
% a  == left endpoint of interval
% b  == right endpoint of interval
% m  == number of subintervals
% f  == function handle for integrand (vectorised)
% OUTPUT:
% Im == computed approximation of integral
h = (b-a)/m;
x = [a+h/2:h:b];
dim = max(size(x));
y = f(x); 
if size(y)==1, y=diag(ones(dim))*y; end;
Im = h*sum(y);

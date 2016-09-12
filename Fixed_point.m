function Fixed_point

g = @(x) x;
% Interval 
a = 1;
b = 2;
% Initial guess
x0 = -1;
% Tolerance
Tol = 1e-8;
% Function Calling
%BMroot = bisect(f,a,b,Tol)
%Fzeroot = fzero(f,x0)
FProot = fpi(g,x0,100)
    

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
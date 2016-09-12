function False_position
f = @(x)  exp(x-1) - x.^3 - 2;
exact = fzero(f,1)
% Function Calling
[br,vecbr] = BM(f,-2,-1,3)

function [xc,x] = BM(f, a, b,Maxit)
if (sign(f(a))*sign(f(b)) > 0);
error('BM condition not satisfied.')
return
end

 Niter_BM = 0;
for k = 1 : (Maxit -1)
c = (b*f(a) - a*f(b))/(f(a) - f(b));
x(k) =c;
%cond = abs(f(c)); % a number to compare with tolerance
if f(c) == 0; xc = c; break; end % solution is found
if sign(f(c))*sign(f(b))< 0; a = c; else b = c; end % get new interval contaning solution
end
xc = c; %new midpoint is best estimate
x(k+1)=xc;
 %Niter_BM = Niter_BM+1
end
end

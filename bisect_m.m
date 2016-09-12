function bisect_m
f = @(x) exp(x-1) - x.^3 - 2;
exact = fzero(f,1)
[br,vecbr] = BM(f,-2,-1,3)
%g = @(x)  (x^3+80)/(2*x^2);
g = @(x)(exp(x)-sin(x))/3;
%[bfp, vecbfp] = FP(f,0,5,1e-5)
%[p1,p2,err,P1] = MM(f,0,.5,1,1e-5,1e-5,1)
%[bsm, vecbsm] = SM(f,0,1,1e-5,5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xc,x] = BM(f, a, b,Maxit)
if (sign(f(a))*sign(f(b)) > 0);
error('BM condition not satisfied.')
return
end
Niter_BM = 0;
for k = 1 : (Maxit -1)
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
end
% Iterative Methods: Jacobi and G-S for large matrix
% Matrix Decomposition: Cholesky factorization
% Objective: Learn how to handle large tridiagonal matrix 
clear; close all; clc
n = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = ones(1,n);
b = [2;v(1:n-2)';2];
D = diag(3*v);
L = -diag(v(1:n-1),-1);
U = -diag(v(1:n-1),1);
% 
Tj = (D)\(-L-U); % Tj for Jacobi
cj = D\b;
% G-S
Tg = -(D+L)\U; % Tj for Gauss-Siedel
cg = (D+L)\b;
x(:,1) = zeros(n,1);
y(:,1) = zeros(n,1);
% Iteration
for k =1:6000
    %Jacobi
    x(:,k+1) = Tj*x(:,k) + cj;
    if norm(x(:,k+1)-v',inf)/norm(v',inf)<1e-6;Jtol=k, break; end
    
end
for k=1:6000
    % G-S
    y(:,k+1) = Tg*y(:,k) + cg;
    if norm(x(:,k+1)-v',inf)/norm(v',inf)<1e-6; Gtol  = k, break; end
end
x=x(:,end)'
y=y(:,end)'
%%%%%%%%%%%%%%%%%%%%% Economized Iterations %%%%%%%%%%%%%%%%%%%%%%%%%
% Jacobi
a = [-1,3,-1];
b = [2,ones(1,n-2),2];
xold = zeros(1,n);
xexact=ones(1,n);

tol =100;
while tol>1e-6
for k = 1:n
    
    if k==1 
        xnew(k) = (-a(3)*xold(k+1) + b(k))/a(2);
    elseif k==n
        xnew(k) =(-a(1)*xold(k-1) + b(k))/a(2);
    else
        xnew(k) = (-a(1)*xold(k-1)-a(3)*xold(k+1) + b(k) )/a(2);
    end
    if k==n;xold = xnew;end
    if k==n; tol =  norm(xexact-xnew,inf)/norm(xexact,inf); end
    
end
end
Jacobi_sol = xnew;
% G - S
xold=0;xnew = xold;
xold = zeros(1,n);
xexact=ones(1,n);
tol =100;
while tol>1e-6
for k = 1:n
    if k==1 
        xnew(k) = (-a(3)*xold(k+1) + b(k))/a(2);
        xold(k) = xnew(k);
    elseif k==n
        xnew(k) =(-a(1)*xold(k-1) + b(k))/a(2);
        xold(k) = xnew(k);
    else
        xnew(k) = (-a(1)*xold(k-1)-a(3)*xold(k+1) + b(k) )/a(2);
        xold(k) = xnew(k);
    end
    if k==n;xold = xnew;end
    if k==n; tol =  norm(xexact-xnew,inf)/norm(xexact,inf); end
    
end
end
G_Sol = xnew;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Chol Factors
%A = R'*R;
B = [5,1,3;1 2 2;3 2 5];
rhs= [4;-1;0];
R=chol(B);
y = R'\rhs; 
x = R\y;
x






















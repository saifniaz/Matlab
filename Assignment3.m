%% Assignment 3
%% Md. Saif Niaz
%% Id: 100555440

format short
A = [1 2 -2;1 1 1; 2 2 1]; b = [7; 2; 5];

%Question 3(a)

D = diag(diag(A));
L = tril(A,-1); 
U = triu(A,1);

Tj = (D)\(-L-U); % Tj for Jacobi
EigTj = eigs(Tj); % Eigenvalues of Tj
SpTj = max(abs(EigTj)) % Spectral radiu
cj = D\b;

%3(C)

Tg = -(D+L)\U; % Tj for Gauss-Siedel
EigTg = eigs(Tg); % Eigenvalues of G-S
SpTg = max(abs(EigTg)) %Spectral radius
cg = (D+L)\b;
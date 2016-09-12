%% Programming Test 2, Question 2
%% Md. Saif Niaz
%% Id: 100555440
clear; close all; clc
format long
n = 3;
A = [3 -1 1; 1 -3 -1; 1 3 5]; b = [0; -4; 2];
e_sol = [1; 2; -1];
v = ones(1,n);
D = diag(diag(A)); % Diagonal Matrix
L = tril(A,-1); % Lower triangular Matrix
U = triu(A,1); % Upper triangular Matrix

Tg = -(D+L)\U % Tj for Gauss-Siedel
EigTg = eigs(Tg) % Eigenvalues of G-S
SpTg = max(abs(EigTg)) %Spectral radius
cg = (D+L)\b;
y(:,1) = [0;0;0];

for k = 1:n
    % G-S
    y(:,k+1) = Tg*y(:,k) + cg;
    if norm(y(:,k+1)-v',inf)/norm(v',inf)<1e-7;Jtol=k, break; end
end
disp('     G - S')
y=y';
sol = [y]

%% Absolute Error with norms 2 is done on paper 

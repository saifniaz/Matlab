% Tutorial5
%Example
% 4x -y +z = 7; 4x -8y +z = -21; -2x + y +5z = 15;
% Do Scaled partial pivoting
% Do Jacobi & G-S and find Spectral radius
% Calc. norm ||.||2 , ||.||inf
format long

A = [3 -1 1; 1 -3 -1; 1 3 5]; b = [0; -4; 2];
    % Jacobi
D = diag(diag(A)); % Diagonal Matrix
L = tril(A,-1); % Lower triangular Matrix
U = triu(A,1); % Upper triangular Matrix

% 
Tj = (D)\(-L-U); % Tj for Jacobi
EigTj = eigs(Tj); % Eigenvalues of Tj
SpTj = max(abs(EigTj)) % Spectral radius
cj = D\b;
% G-S
Tg = -(D+L)\U; % Tj for Gauss-Siedel
EigTg = eigs(Tg); % Eigenvalues of G-S
SpTg = max(abs(EigTg)) %Spectral radius
cg = (D+L)\b;
x(:,1) = [0;0;0];
y(:,1) = [0;0;0];
% Iteration
for k = 1:10
    %Jacobi
    x(:,k+1) = Tj*x(:,k) + cj;
    % G-S
    y(:,k+1) = Tg*y(:,k) + cg;
end
x=x';
y=y';
disp('Jacobi                               G-S')
sol = [x,y]

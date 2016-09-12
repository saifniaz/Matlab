% Tutorial5
%Example
% 4x -y +z = 7; 4x -8y +z = -21; -2x + y +5z = 15;
% Do Scaled partial pivoting
% Do Jacobi & G-S and find Spectral radius
% Calc. norm ||.||2 , ||.||inf

n = 10;
v = ones(1, n);
b = [2;v(1:n-2)';2];
D = diag(3*v);
L = -diag(v(1:n-1), -1);
U = -diag(v(1:n-1), 1);

%A = [4 -1 1;4 -8 1; -2 1 5]; b = [7;-21;15];
    % Jacobi
%D = diag(diag(A)); % Diagonal Matrix
%L = tril(A,-1); % Lower triangular Matrix
%U = triu(A,1); % Upper triangular Matrix

% 
Tj = (D)\(-L-U); % Tj for Jacobi
%EigTj = eigs(Tj); % Eigenvalues of Tj
%SpTj = max(abs(EigTj)) % Spectral radius
cj = D\b;
% G-S
Tg = -(D+L)\U; % Tj for Gauss-Siedel
%EigTg = eigs(Tg); % Eigenvalues of G-S
%SpTg = max(abs(EigTg)) %Spectral radius
cg = (D+L)\b;
x(:,1) = zeros(n,1);
y(:,1) = zeros(n,1);
% Iteration
for k =1:6000
    %Jacobi
    x(:,k+1) = Tj*x(:,k) + cj;
    if norm(x(:,k+1)-v',inf)/norm(v',inf)<1e-6;Jtol=k, break; end
    % G-S
    %y(:,k+1) = Tg*y(:,k) + cg;
end
for k =1:6000
    %Jacobi
    y(:,k+1) = Tg*y(:,k) + cg;
    if norm(x(:,k+1)-v',inf)/norm(v',inf)<1e-6;Jtol=k, break; end
    % G-S
    %y(:,k+1) = Tg*y(:,k) + cg;
end
x=x(:,end)'
y=y(:,end)'
% disp('Jacobi                               G-S')
% sol = [x,y]

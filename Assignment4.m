%% Assignment 4
%% Md. Saif Niaz
%% Id: 100555440

p = @(x) 0.97326 + 3.8438 * x - 4.905*x.^2;

A = [0.97326 0.0265 2.3329e-04; 0.97326 3.8338 4.8795; 0.97326 0.0146 7.0828e-05; 0.97326 3.8438 4.905]

b = [0.9995; 3.8338; 0.0146; -0.0879]

[Q,R] = qr(A)
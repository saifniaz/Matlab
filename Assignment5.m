%% Assignment 5
%% Md.Saif Niaz
%% 100555440

%Question 2(a)--------------------------------------------------
y1 = [3 ; 1; 2];
v = [1 -1 1;
     1 0 0;
     1 1 1];
  
a = v\y1

%Question 2 (b)-------------------------------------------------

x = [ -1.0, 0.0, 1.0 ];
y = [ 3.0, 1.0, 2.0 ];
Pi2 = polyfit( x, y, 2 ) 
xx = linspace( -2, 2 ); 
yy = polyval( Pi2, xx); 
plot( x, y, 'ro', xx, yy, 'b-' ) 
grid on
axis([-5 5 0 7])

 
 
function simpplot(M)

% function to illustrate Composite Simpson's rule and error graphically

% close all
f=@humps;
a=0;
b=1;
numpoints=60*M+1;
%if M<3, numpoints=61, end
x = linspace(a,b,numpoints);
y = f(x);
xx = [ a x b ];
yy = [ 0 y 0 ];

figure(2)
clf
title('Error in Composite Simpson''s Rule')
figure(1)
clf
title('Composite Simpson''s Rule')
line(x,y,'color','b')
xl = a;
for k=1:M
    xr = xl + (b-a)/M;
    ind = find( (x<=xr) & (x>=xl) );
    xx = [xl x(ind) xr];
    xq = [ xl xl+0.5*(b-a)/M xr ];
    yq = f(xq);
    p = polyinterp( xq, yq, x(ind) );
    pp = [ 0 p 0 ];
    figure(1)
    hp = patch(xx,pp,'g');
    set(hp,'FaceAlpha',0.2)
    Qpts = line(xq,yq);
    set(Qpts,'Color','b','LineStyle','none','MarkerFaceColor','b','Marker','o')
    figure(2)
    xx = [ x(ind) fliplr(x(ind)) ];
    ee = [ y(ind) fliplr(p) ];
    hp = patch(xx,ee,'c');
    set(hp,'FaceAlpha',0.2)
    xl=xr;
end

function yy = polyinterp(x,y,xx)
w = barywts(x);
yy = baryval(x,y,w,xx);
end

function w = barywts( x )
% w = barywts( x )
% D. Aruliah, Oct. 2008
% Function to compute barycentric weights associated with interpolation
% nodes x.
% INPUT:
% x -- distinct interpolation nodes, vector of length n
% OUTPUT
% w -- barycentric weights based on nodes x, vector of length n
% See Berrut & Trefethen (SIREV, 2004)
% Note: At present, contains no verification to ensure distinct nodes

% Calculation is vectorised to speed up loops
n = length( x );
w = repmat(shiftdim(x), 1, n );
w = prod( w' - w + eye(n) ); % Do all products/differences
w = 1./w;
end

function yy = baryval( x, y, w, xx )
% yy = baryval( x, y, w, xx )
% D. Aruliah, Oct. 2008
% INPUT:
% x  -- distinct interpolation nodes, vector of length n
% y  -- coordinates, vector of length n
% w  -- barycentric weights (precomputed), vector of length n
% xx -- array of values on which to evaluate interpolant
% OUTPUT:
% yy -- values of polynomial interpolant on grid defined by xx
% Evaluate polynomial interpolant using barycentric Lagrange interpolation
% (of the first kind)
% See Berrut & Trefethen (SIREV, 2004)
% Note: At present, contains no verification to ensure distinct nodes

% Calculation is vectorised to speed up loops
n = length(y);
yy = zeros(size(xx));
omega = ones(size(xx));        % array to store values of nodal polynomial
exact = yy;
for k=1:n
    diffx = xx-x(k);
    exact( diffx==0 ) = k;               % Trick from Berrut & Trefethen
    yy = yy + w(k) * y(k) ./ diffx;
    omega = omega .* diffx;
end
yy = yy .* omega;
ii = find( exact );    % Determine which points coincide with nodes
yy(ii) = y(exact(ii)); % Substitute in values from interpolation data
end

end

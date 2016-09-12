% interp_test.m
% Routine to plot and quantify convergence of piecewise polynomial
% approximations to a prescribed function

% Define function to approximate and finest grid
f = @(x) sin(exp(3*x));
xx = linspace( -1, 1, 5000 );
yy = f( xx );

disp( [ 'Available methods: (n) nearest, (l) linear, ' ...
        '(p) pchip, (s) spline' ])
method = input('Please enter a choice of interpolation methods: ', 's');

% Initialise empty array to accumulate results.
Results = [];
figure(1), clf
set(gca,'FontSize',16)
hf = line(xx,yy);
set(hf,'Color','k','LineStyle','--')
hd = line( 0, 0, 'LineStyle','none','Marker','x');
hS = line( 0, 0, 'Color','r','LineWidth',2,'LineStyle','-.' );
text(-0.75,-0.5,'$f(x)=\sin(e^{3x})$','Interpreter','latex','FontSize',16);
axis([-1 1 -1 1]), axis equal, axis tight
% Start looping
for m = 2.^[ 1:11 ]
    h=1/m;
    t = linspace(-1,1,m+1); % Define spline knots
    y = f(t); % Sample at knots
    switch method
        case {'n','N'}
            t_tmp = t;
            y_tmp = y;
            t = 0.5*( t(1:end-1) + t(2:end) );
            t = [ -1-h, t, 1+h ];
            y = f(t);
            S = interp1( t, y, xx, 'nearest' );
            t = t_tmp;
            y = y_tmp;
        case {'l','L'}
            S = interp1( t, y, xx, 'linear' );
        case {'p','P'}
            S = pchip( t, y, xx );
        case {'s','S'}
            S = spline( t, y, xx );
    end
    
    if strcmp(method,'nearest'),  end
    err = norm( yy-S , 'inf');
    set(hd,'Xdata',t,'Ydata',0*t)
    set(hS,'Xdata',xx,'Ydata', S);
    ht = title(['$m=' num2str(m) '\ ,\|{f}{{-}}{S}\|_\infty =$ ' ...
        num2str(err,'%10.3e')]);
    set(ht,'Fontsize',16,'Interpreter','latex')
    Results = [ Results, [h; err; 0] ];
    disp('Press return to continue.');
    pause
end
Results(3,2:end) = log( Results(2,1:end-1)./Results(2,2:end) )/log(2);
S1 = sprintf('    h\t\t Err\t\t     p\n');
disp( [S1, sprintf('%6.2e\t%10.4e\t%8.2f\n', Results) ] )

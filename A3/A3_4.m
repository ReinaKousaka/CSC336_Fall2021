x = [ 0.9  1.3  1.9  2.1  2.6  3.0  3.9  4.4  4.7  5.0  6.0  7.0  8.0  9.2 ...
     10.5 11.3 11.6 12.0 12.6 13.0 13.3];
y = [ 1.3  1.5  1.85 2.1  2.6  2.7  2.4  2.15 2.05 2.1  2.25 2.3  2.25 1.95 ...
      1.4  0.9  0.7  0.6  0.5  0.4  0.25];
xi  = linspace(0.87, 13.33, 100);

% not-a-know cubic spline interpolant
yvs = spline(x, y, xi);

% 20th degree polynomial interpolant
p = polyfit(x, y, 20);
polyfits = polyval(p, xi);

% plot
plot(xi, yvs, 'k-', ...
    xi, polyfits, '--', ...
    x, y, 'ro');
axis([0 14 -6 5]);
legend('Cubic Spline', 'Polynomial');
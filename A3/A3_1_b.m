% computing the spreading of influenza
% test for Newton's on one (large) timestep

% beta transmission, gamma recovery, mu death/birth (replenishment)
beta = 0.75; gamma = 0.06; mu = 0.01/365;
% rho vaccination, 1/omega, 1/omegaV immunity from recovering, vaccination
rho = 150 * mu; omega = 2/365; omegaV = 1.50/365;

% initial conditions
N = 15e6;
h = 1;
y02 = 3800; y03 = 589533; y04 = 0.80*N;
y01 = N - y02 - y03 - y04;
y0 = [y01 y02 y03 y04]';

dt = 1; % stepsize for time is h (or dt)
Beta = dt*beta/N; Gamma = dt*gamma; Mu = dt*mu; % for convenience
Rho = dt*rho; Omega = dt*omega; OmegaV = dt*omegaV;

maxit = 10; tol = 1e-7; % Newton's parameters
fprintf(' k   S         I       R        V        Total     Residual\n');
y = y0;
yinit = y; % initial guess for Newton's
for k = 1:maxit
    % define vector f and its inf norm
%     f = [
%         y(1) - yinit(1) - h * mu * N + h * mu * y(1) + h * beta * y(1) * y(2) / N + h * rho * y(1) - h * omega * y(3) - h * omegaV * y(4);
%         y(2) - yinit(2) + h * mu * y(2) - h * beta * y(1) * y(2) / N + h * gamma * y(2);
%         y(3) - yinit(3) + h * mu * y(3) - h * gamma * y(2) + h * omega * y(3);
%         y(4) - yinit(4) + h * mu * y(4) - h * rho * y(1) + h * omegaV * y(4);
%     ];
    f = [
        y(1) - yinit(1) - Mu * N + Mu * y(1) + Beta * y(1) * y(2) + Rho * y(1) - Omega * y(3) - OmegaV * y(4);
        y(2) - yinit(2) + Mu * y(2) - Beta * y(1) * y(2) + Gamma * y(2);
        y(3) - yinit(3) + Mu * y(3) - Gamma * y(2) + Omega * y(3);
        y(4) - yinit(4) + Mu * y(4) - Rho * y(1) + OmegaV * y(4);
    ];
    fnorm = norm(f, inf);
    fprintf('%2d %9.0f %6.0f %9.0f %6.0f %10.0f %9.2e\n', ...
        k-1, y, sum(y), fnorm);
    % stopping criterion
    if fnorm <= tol
        break;
    end
    % define Jacobian matrix
    J = [
        Mu + Rho + 1 + Beta * y(2), Beta * y(1), -Omega, -OmegaV;
        -Beta * y(2), Mu + Gamma + 1 - (Beta * y(1)), 0, 0;
        0, -Gamma, Mu + Omega + 1, 0;
        -Rho, 0, 0, Mu + OmegaV + 1
    ];
    s = J \ f;
%         fprintf('s is %f %f %f %f', s);
    % apply Newton's iteration to compute new y
    y = y - s;
end
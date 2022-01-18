% beta transmission, gamma recovery, mu death/birth (replenishment)
beta = 0.25; gamma = 0.06; mu = 0.01/365;
%1/omega, 1/omegaV immunity from recovering, vaccination
omega = 2/365; omegaV = 1.50/365;
% initial conditions
N = 15e6;
h = 1 / 24;
y02 = 3800; y03 = 589533; y04 = 0.80*N;
y01 = N - y02 - y03 - y04;
y0 = [y01 y02 y03 y04]';

dt = 1/24; % stepsize for time is h (or dt)
Beta = dt*beta/N; Gamma = dt*gamma; Mu = dt*mu; % for convenience
Omega = dt*omega; OmegaV = dt*omegaV;

t_end = 4 * 365 + 1;
nstep = t_end / h;

% difference   
num_Newton = zeros(4, nstep); % store the number of Newton's iteration
for rho_index = 1:4
    rho = rho_index * 150 * mu;
    Rho = dt * rho;
    fprintf('rho = %d * mu\n', rho_index * 150);
    % result matrix
    yi = zeros(nstep + 1, 4);   % row i <---> time at i-1
    yi(1, :) = y0;
    fprintf(' k      S         I       R        V        Total\n');
    fprintf('%d %9.0f %6.0f %9.0f %6.0f %10.0f   Initial\n', rho_index * 150, y0, sum(y0));
    for i = 1:nstep
        maxit = 10; tol = 1e-7; % Newton's parameters
    
        y = yi(i, :);  % the initial guess is the solution of the previous time
        for k = 1:maxit       
            % define vector f and its inf norm
               f = [
                y(1) - yi(i, 1) - Mu * N + Mu * y(1) + Beta * y(1) * y(2) + Rho * y(1) - Omega * y(3) - OmegaV * y(4);
                y(2) - yi(i, 2) + Mu * y(2) - Beta * y(1) * y(2) + Gamma * y(2);
                y(3) - yi(i, 3) + Mu * y(3) - Gamma * y(2) + Omega * y(3);
                y(4) - yi(i, 4) + Mu * y(4) - Rho * y(1) + OmegaV * y(4);
            ];
            fnorm = norm(f, Inf);
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
            s = J \ (-f);
            % apply Newton's iteration to compute new y
            y = y + s.';
        end
        yi(i + 1, :) = y;
        num_Newton(rho_index, i) = k - 1;
    end
    % done for one rho, output
    fprintf('%d %9.0f %6.0f %9.0f %6.0f %10.0f   Final\n', rho_index * 150, yi(nstep + 1, :), sum(yi(nstep + 1, :)));
    
    fprintf('%d %9.0f %6.0f %9.0f %6.0f    MAX\n', rho_index * 150, max(yi(1:nstep + 1, 1)), max(yi(1:nstep + 1, 2)),max(yi(1:nstep + 1, 3)),max(yi(1:nstep + 1, 4)));
    [argvalue, argmax] = max(yi(:, 2));
    fprintf('The max # infected occurs on day %d\n', (argmax - 1)*h);
    fprintf('\n');
    % plot
    X = ((1:nstep + 1) - 1) * h;
    figure();
    plot(X, yi(:, 1), ':k', ...
        X, yi(:, 2), '--r', ...
        X, yi(:, 3), '-.b', ...
        X, yi(:, 4), '-g');
    legend('Susceptible', 'Infected', 'Recovered', 'Vaccinated');
    xlabel('Time Point');
    ylabel('Number of People'); 
    axis tight;
end 
% plot the 5th figure
figure();
plot((1:nstep), num_Newton(1, :), ':k', ...
     (1:nstep), num_Newton(2, :), '--r', ...
     (1:nstep), num_Newton(3, :), '-.b',...
     (1:nstep), num_Newton(4, :), '-g');
legend('150', '300', '450', '600');
xlabel('Time Point');
ylabel('Number of Newton Iteration');
axis([1 nstep 0.9 2.1]);

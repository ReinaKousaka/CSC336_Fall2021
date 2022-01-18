index = 0;
ni = [4 8 16 32];
for n = ni
    N = 3*n-2;
    A = speye(N, N);
    % set loop equations
    A(1, 1:3) = [1 1 1];
    A(N, N-2:N) = [-1 0 2];
    for i = 4:3:(3*n-3)
        A(i, i-2:i+2) = [-1 0 1 1 1];
    end
    % set node equations
    for i = 1:3:(3*n-5)
        A(i+1, i:i+3) = [1 -1 0 -1];
        if i ~= 3*n-5
            A(i+2, i+1:i+5) = [1 -1 0 0 1];
        else
            A(N-1, N-2:N) = [1 -1 1]; 
        end
    end
    % generate RHS vector
    b = zeros(N, 1);
    b(1) = 100;
    % solve the system
    x = A\b;
    % store info for plotting
    index = index + 1;
    for i = 1:n
        t(i, index) = x(i * 3 - 2);
    end
    % output
    fprintf('n = %d\n', n);
    fprintf('Max Intensity: %.3f\n', max(x));
    fprintf('Min Intensity: %.3f\n', min(x));
    fprintf('Sum of Intensities: %.3f\n', sum(x));
    fprintf('Condition Number: %.3f\n', condest(A));
    if n == 8 
        [L, U] = lu(A);
        figure(1);
        spy(A);
        figure(2);
        spy(L);
        figure(3);
        spy(U);
    end
end
figure(4);
plot([1:ni(1)]/ni(1), t(1:ni(1), 1), 'r-', ...
[1:ni(2)]/ni(2), t(1:ni(2), 2), 'g--', ...
[1:ni(3)]/ni(3), t(1:ni(3), 3), 'b-.', ...
[1:ni(4)]/ni(4), t(1:ni(4), 4), 'k.');
legend('n=4', 'n=8', 'n=16', 'n=32')
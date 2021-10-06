% ln'(1) = 1,  ln''(1) = -1
h = 10.^ (-16:-1);
x = ones(1, 16);

g_a = (log(x + h) - log(x)) ./ h;
g_b = (log(x + h) - log(x - h)) ./ (2 * h);
g_c = (log(x + h) + log(x - h) - 2 * log(x)) ./ (h.^2);   
err_a = (1 - g_a) / 1;
err_b = (1 - g_b) / 1;
err_c = (-1 - g_c) / (-1);

epsilon = eps .* ones(1, 16);
m = ones(1, 16);
bound_a = (m .* h) / 2 + 5 * epsilon ./ h;
bound_b = (2 * m) .* (h.^2) / 6 + 5 * epsilon ./ (2 * h);
bound_c = (6 * m) .* (h.^2) / 12 + 6 * epsilon ./ (h.^2);

loglog(h, abs(err_a), '-', h, bound_a, 'x-', ...
    h, abs(err_b), '--', h, bound_b, 'x--', ...
    h, abs(err_c), ':', h, bound_c, 'x:');
legend('err-a', 'bound-a', 'err-b', 'bound-b', 'err-c', 'bound-c');
axis tight;

[~, index_a] = min(abs(err_a));
[~, index_b] = min(abs(err_b));
[~, index_c] = min(abs(err_c));

fprintf('Min Error in case (a): stepsize = %9.2e, error = %9.2e\n', ...
    h(index_a), err_a(index_a));
fprintf('Min Error in case (b): stepsize = %9.2e, error = %9.2e\n', ...
    h(index_b), err_b(index_b));
fprintf('Min Error in case (c): stepsize = %9.2e, error = %9.2e\n', ...
    h(index_c), err_c(index_c));
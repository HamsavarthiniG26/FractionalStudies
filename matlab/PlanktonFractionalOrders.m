Lambda = 2.8;
beta = 0.07;
a0 = 0.05;
a1 = 0.05;
a2 = 0.06;
a3 = 0.3;
B = 0.15;
alpha = 0.1;
l = 0.07;
g = 0.99;
zeta = 0.01;
delta = 0.07;
t0 = 0;
tf = 50;
h = 0.01;
t = t0:h:tf;
C0 = 0.99;
U0 = 0.9;
V0 = 0.5;
orders = [0.6, 0.7, 0.8, 0.9];
results_C = zeros(length(orders), length(t));
results_U = zeros(length(orders), length(t));
results_V = zeros(length(orders), length(t));
results_C(:, 1) = C0;
results_U(:, 1) = U0;
results_V(:, 1) = V0;
for idx = 1:length(orders)
 q = orders(idx);
 N = length(t);
 w = zeros(1, N);
 for k = 0:N-2
 arg1 = q + 1;
 arg2 = k + 1;
 arg3 = q - k + 1;
 if arg3 > 0
 w(k + 1) = (-1)^k * gamma(arg1) / (gamma(arg2) * gamma(arg3));
 else
 w(k + 1) = 0;
 end
 end
 for n = 2:N
 C_sum = 0; U_sum = 0; V_sum = 0;
 for k = 0:n-2
 if n - 1 - k > 0 && (q - (n - 1 - k) + 1) > 0
 C_sum = C_sum + w(k + 1) * (results_C(idx, n - 1 - k) - 
results_C(idx, 1));
 U_sum = U_sum + w(k + 1) * (results_U(idx, n - 1 - k) - 
results_U(idx, 1));
 V_sum = V_sum + w(k + 1) * (results_V(idx, n - 1 - k) - 
results_V(idx, 1));
 end
 end
 dCdt = @(C, U, V) Lambda * (1 - (C ./ (C + a0))) .* U - (U .* C ./ (C 
+ a2)) - beta .* C .* V ./ (C + a3) - C;
 dUdt = @(C, U, V) (B .* C ./ (C + a1) - alpha .* U) .* U - (g .* U .* 
V ./ (U + l)) - delta .* U;
 dVdt = @(C, U, V) (beta .* C.^2 ./ (C.^2 + a2^4)) .* U .* V ./ (U + 
l) - zeta .* V;
 results_C(idx, n) = results_C(idx, n - 1) + h^q * dCdt(results_C(idx, 
n - 1), results_U(idx, n - 1), results_V(idx, n - 1)) + h^q * C_sum;
 results_U(idx, n) = results_U(idx, n - 1) + h^q * dUdt(results_C(idx, 
n - 1), results_U(idx, n - 1), results_V(idx, n - 1)) + h^q * U_sum;
 results_V(idx, n) = results_V(idx, n - 1) + h^q * dVdt(results_C(idx, 
n - 1), results_U(idx, n - 1), results_V(idx, n - 1)) + h^q * V_sum;
 end
end
figure;
hold on;
for idx = 1:length(orders)
 plot(t, results_C(idx, :), 'LineWidth', 1.5, 'DisplayName', 
sprintf('Order=%.1f', orders(idx)));
end
xlim([0 10]);
xlabel('Time');
ylabel('C(t)');
legend('show');
hold off;
figure;
hold on;
for idx = 1:length(orders)
 plot(t, results_U(idx, :), 'LineWidth', 1.5, 'DisplayName', 
sprintf('Order=%.1f', orders(idx)));
end
xlim([0 10]);
xlabel('Time');
ylabel('U(t)');
legend('show');
hold off;
figure;
hold on;
for idx = 1:length(orders)
 plot(t, results_V(idx, :), 'LineWidth', 1.5, 'DisplayName', 
sprintf('Order=%.1f', orders(idx)));
end
xlabel('Time');
ylabel('V(t)');
legend('show');
hold off;

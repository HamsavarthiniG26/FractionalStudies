% Parameters
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
h = 0.1;
t = t0:h:tf;
C0 = 0.99;
U0 = 0.9;
V0 = 0.5;
qs = [0.6, 0.7, 0.8, 0.9]; % Different values of q

% Loop through each order value of q
for q_idx = 1:length(qs)
    q = qs(q_idx);
    
    % Initialize results
    results_C = zeros(1, length(t));
    results_U = zeros(1, length(t));
    results_V = zeros(1, length(t));
    results_C(1) = C0;
    results_U(1) = U0;
    results_V(1) = V0;
    N = length(t);
    w = zeros(1, N);
    
    % Calculate fractional weights
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
    
    % Time stepping
    for n = 2:N
        C_sum = 0; U_sum = 0; V_sum = 0;
        for k = 0:n-2
            if n - 1 - k > 0 && (q - (n - 1 - k) + 1) > 0
                C_sum = C_sum + w(k + 1) * (results_C(n - 1 - k) - results_C(1));
                U_sum = U_sum + w(k + 1) * (results_U(n - 1 - k) - results_U(1));
                V_sum = V_sum + w(k + 1) * (results_V(n - 1 - k) - results_V(1));
            end
        end
        % Differential equations
        dCdt = @(C, U, V) Lambda * (1 - (C / (C + a0))) * U - (U * C / (C + a2)) ...
                          - beta * C * V / (C + a3) - C;
        dUdt = @(C, U, V) (B * C / (C + a1) - alpha * U) * U - (g * U * V / (U + l)) - delta * U;
        dVdt = @(C, U, V) (beta * C^2 / (C^2 + a2^4)) * U * V / (U + l) - zeta * V;
        
        % Update results using fractional method
        results_C(n) = results_C(n - 1) + h^q * dCdt(results_C(n - 1), ...
                          results_U(n - 1), results_V(n - 1)) + h^q * C_sum;
        results_U(n) = results_U(n - 1) + h^q * dUdt(results_C(n - 1), ...
                          results_U(n - 1), results_V(n - 1)) + h^q * U_sum;
        results_V(n) = results_V(n - 1) + h^q * dVdt(results_C(n - 1), ...
                          results_U(n - 1), results_V(n - 1)) + h^q * V_sum;
    end
    
    % Plot results for this value of q
    figure;
    plot3(results_C, results_U, results_V, 'LineWidth', 1.5);
    xlabel('C(t)');
    ylabel('U(t)');
    zlabel('V(t)');
    title(['3D Plot for q = ', num2str(q)]);
    grid on;
    view(3);
end

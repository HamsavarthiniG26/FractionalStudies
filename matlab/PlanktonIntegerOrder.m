Lambda = 2.8; % Value of Lambda
beta = 0.07; % Value of beta
a0 = 0.05; % Value of a0
a1 = 0.05; % Value of a1
a2 = 0.06; % Value of a2
a3 = 0.3; % Value of a3
B = 0.15; % Value of B
alpha = 0.1; % Value of alpha
l = 0.07; % Value of l
gamma = 0.99; % Value of gamma
zeta = 0.01; % Value of zeta
delta = 0.07; % Value of delta
% Initial conditions
C0 = 0.99; % Initial condition for C
U0 = 0.9; % Initial condition for U
V0 = 0.5; % Initial condition for V
initial_conditions = [C0; U0; V0];
% Time span for simulation
tspan = [0 10];
[t, X] = ode45(@(t, X) system(t, X, Lambda, beta, a0, a1, a2, a3, B, alpha, 
l, gamma, zeta, delta), tspan, initial_conditions);
figure;
plot(t, X(:, 1), 'r', 'LineWidth', 2); % Plot for C
hold on;
plot(t, X(:, 2), 'g', 'LineWidth', 2); % Plot for U
plot(t, X(:, 3), 'b', 'LineWidth', 2); % Plot for V
hold off;
ylabel('Population/Dynamics');
legend('C(t)', 'U(t)', 'V(t)', 'Location', 'best');
grid on;
function dXdt = system(~, X, Lambda, beta, a0, a1, a2, a3, B, alpha, l, 
gamma, zeta, delta)
 C = X(1);
 U = X(2);
 V = X(3);
 dCdt = Lambda * (1 - (C / (C + a0))) * U - (U * C / (C + a2)) - beta * C 
* V / (C + a3) - C;
 dUdt = (B * C / (C + a1) - alpha * U) * U - (gamma * U * V / (U + l)) - 
delta * U;
 dVdt = (beta * C^2 / (C^2 + a2^4)) * U * V / (U + l) - zeta * V;
 dXdt = [dCdt; dUdt; dVdt];
end
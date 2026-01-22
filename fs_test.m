clear; clc;

% -------------------------
% PARÁMETROS DEL MATERIAL
% -------------------------
phic_deg = 32;       % Ángulo de fricción crítica [grados]
phic = deg2rad(phic_deg); % Convertir a radianes

eta = 0.0;           % Parámetro de forma (0 = Drucker-Prager)
a01 = real(1 / exp(0)); % a0*eta = 0 => exp(0) = 1 => a01 = 1
ag = sqrt(3)/2 * (3 - sin(phic)) / (sqrt(2)*sin(phic)) * a01;

% Presión media constante
p = 100; % [kPa]

% -------------------------
% CÁLCULO DE LA ENVOLVENTE
% -------------------------
N = 100;
theta = linspace(-pi/6, pi/6, N);  % Ángulos de Lode
p_vals = zeros(N,1);
q_vals = zeros(N,1);

M = [0 1 1; 1 0 1; 1 1 0];  % Para calcular I2

for i = 1:N
    t = theta(i);

    % Tensiones desviadoras normalizadas
    s1 = sqrt(2/3) * cos(t);
    s2 = sqrt(2/3) * cos(t - 2*pi/3);
    s3 = sqrt(2/3) * cos(t + 2*pi/3);
    s_vec = [s1; s2; s3];

    % Tensor de tensiones principales
    sigma_princ = p * (1 + ag * s_vec);

    % Cálculo de invariantes
    I1 = sum(sigma_princ);
    I2 = 0.5 * (sigma_princ' * M * sigma_princ - I1^2);
    I3 = sigma_princ(1)*sigma_princ(2)*sigma_princ(3);

    q = sqrt(1.5 * (sigma_princ' * M * sigma_princ - I1^2 / 3));

    % Almacenar resultados
    p_vals(i) = p;
    q_vals(i) = q;
end

% -------------------------
% PLOT
% -------------------------
figure;
plot(p_vals, q_vals, 'b-', 'LineWidth', 1.5);
xlabel('p [kPa]');
ylabel('q [kPa]');
title('Figura 3 - Envolvente p–q (basado en la Figura 3 de Qian et al., 2024)');
grid on;
axis equal;

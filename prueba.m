clear all
clc;

phi = [linspace(0,1,100)]';
beta = [linspace(0,1,100)]';
epsacc = [linspace(0,1,100)]';

phimod = phi.*(1-beta)+beta.*phi.*epsacc;

plot(phimod,epsacc);

%%
% Parámetros iniciales para la superficie de fluencia en 3D
Mf = 1.68; % Coeficiente de fricción
n = 0.5; % Índice de efecto de presión hidrostática
pr = 100; % Presión de referencia (kPa)
phi_deg = 34; % Ángulo de fricción interna (grados)
phi_rad = deg2rad(phi_deg); % Conversión a radianes
c = 0; % Cohesión inicial (kPa)
a = 0.5; % Parámetro de interpolación entre SMP y Mises

% Definición de la ecuación de la superficie de fluencia según Yao et al. (2015)
% Basado en la ecuación (12) y la ecuación (16) de Yao et al. (2015)
alpha_c_sq = (4 * sin(phi_rad)^2) / (3 * (3 - sin(phi_rad))^2);
B2 = sqrt(alpha_c_sq / 2);

% Rango de tensiones principales con mayor resolución
theta = linspace(0, 2*pi, 100); % Más puntos para suavizar la superficie
p_range = linspace(0, 1000, 100); % Expansión del dominio de p

% Inicialización de matrices
[P, Theta] = meshgrid(p_range, theta);
Q = Mf * (P / pr).^n * pr; % Aplicación de la ecuación (12) de Yao et al.

% Cálculo de los invariantes de tensión correctamente
sigma1 = P + Q .* cos(Theta);
sigma3 = P - Q .* cos(Theta);
sigma2 = P - Q .* sin(Theta);

I1_t = sigma1 + sigma2 + sigma3;
I2_t = sigma1 .* sigma2 + sigma2 .* sigma3 + sigma3 .* sigma1;
I3_t = sigma1 .* sigma2 .* sigma3;

% Interpolación entre SMP y Mises con la ecuación (16)
q_star = a * sqrt(I1_t.^2 - 3*I2_t) + (1-a) * (2/sqrt(3)) * I1_t;

% Graficar la superficie de fluencia en 3D con mejor resolución
figure;
surf(sigma1, sigma3, sigma2, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
xlabel('\sigma_1 (kPa)', 'FontSize', 14);
ylabel('\sigma_3 (kPa)', 'FontSize', 14);
zlabel('\sigma_2 (kPa)', 'FontSize', 14);
title('Superficie de Fluencia en el Espacio de Tensiones (USC)', 'FontSize', 16);
colorbar;
grid on;
view(135,30);
axis equal;

% Referencias:
% - Ecuación (12): Definición de la superficie de fluencia 【74:0†Yao_etal_2015.pdf】
% - Ecuación (14): Cálculo de q basado en los invariantes de esfuerzo 【74:1†Yao_etal_2015.pdf】
% - Ecuación (16): Interpolación entre SMP y Mises 【74:1†Yao_etal_2015.pdf】
clear all
clc

 input_data_biobio;
 init_state_biobio;
 
[SS,EE,INV_S,INV_E,HARD] = updateModel(y0,parms,nspb,path_info);

T = [SS(1), SS(4), SS(5);
     SS(4), SS(2), SS(6);
     SS(5), SS(6), SS(3)];

I1 = trace(T);
I2 = T(1,1)*T(2,2) + T(2,2)*T(3,3) + T(3,3)*T(1,1);
I3 = T(1,1)*T(2,2)*T(3,3);

phi_c = 34; 
alpha_c = sqrt((4 * sind(phi_c)^2) / (3 * (3 - sind(phi_c))^2)); % Eq.(15)
B2 = 1 / (sqrt(2) * alpha_c); % Relación entre B2 y phi_c (Ecuación 16)
B1 = 1 / (2*alpha_c^2);

m = [1;1;1;0;0;0]';
T_norm = T/trace(T);
T_dev = T_norm-eye(3)/3;

gc = 0.5*trace(T_dev^2)-(0.5)*(1/B1^2);
fc = 0.5*trace(T_dev^2)-(0.5)*(1/B2^2);

q_D = sqrt(I1^2-3*I2)
q_M = 2*I1 / (3*sqrt((I1*I2-I3)/(I1*I2-9*I3))-1)
%%

% Datos iniciales y parámetros
T = [100 30 20; 30 80 25; 20 25 60]; % Tensor de tensiones (MPa)
trace_T = trace(T);
T_norm = T / trace_T;
T_dev = T_norm - eye(3)/3; % Parte desviadora del tensor de tensiones

phi_c = 34; % Ángulo de fricción crítica (grados)
alpha_c = sqrt((4 * sind(phi_c)^2) / (3 * (3 - sind(phi_c))^2)); % Eq. (15)
B1 = 1 / (2 * alpha_c^2);
B2 = 1 / (sqrt(2) * alpha_c);

% Tensor de deformaciones iniciales EE
EE = [0.001 0 0; 0 0.0008 0; 0 0 0.0005];

% Cálculo de \mathcal{L} y \mathcal{N}
L = @(T, D) D + trace(D) * T + B1 * trace(T * D) * T; % Ecuación (5)
N = @(T) B2 * (T + T_dev);                        % Ecuación (6)

% Evaluación de \mathcal{L} y \mathcal{N}
L_val = L(T_norm, EE);
N_val = N(T_dev);

% Inversa de \mathcal{L} (simplificación para prueba)
L_inv = inv(L_val); % Nota: Revisar singularidades para casos reales

% Cálculo de \vec{D}_c
D_c = -L_inv * N_val;

% Cálculo de invariantes de esfuerzos
I1 = trace(T);
J2 = 0.5 * trace(T_dev^2);
q = sqrt(3 * J2);
p = I1 / 3;

% Graficar la superficie de fluencia (Ecuación 14)
p_vals = linspace(0, max(p)*1.5, 100);
q_vals = sqrt(2) * alpha_c * p_vals;

figure;
plot(p_vals, q_vals, 'r', 'LineWidth', 2);
hold on;
xlabel('Presión media p [MPa]');
ylabel('Esfuerzo desviador q [MPa]');
title('Superficie de fluencia en el espacio (p, q)');
grid on;
hold off;

% Resultados
disp('Tensor \mathcal{L}:');
disp(L_val);
disp('Tensor \mathcal{N}:');
disp(N_val);
disp('Tensor de deformación en estado crítico (\vec{D}_c):');
disp(D_c);
%%

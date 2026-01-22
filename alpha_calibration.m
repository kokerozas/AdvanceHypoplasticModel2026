%% Script para calcular el parámetro alpha para múltiples ensayos triaxiales
clear all;clc

% -------------------------------
% 1. Importar datos desde Excel
% -------------------------------
filename = 'data_alpha.xlsm'; % Nombre del archivo Excel
sheetName = 'Hoja1'; % Nombre de la hoja
data = readmatrix(filename, 'Sheet', sheetName);

% Índices de las columnas para cada ensayo
% TMD DR VARIATION
ensayo_indices = [1, 2, 3, 4;
                  6, 7, 8, 9;
                 11, 12, 13, 14; 
                 16, 17, 18, 19];
% % TMD S3 VARIATION
% ensayo_indices = [1, 2, 3, 4;
%                   6, 7, 8, 9;
%                  11, 12, 13, 14];

             
num_ensayos = size(ensayo_indices, 1);

% %Cuitiño,2016
e_d = [0.508]; 
e_i = [.947];
e_c = [0.824]; 

% %Fbrown,2024
% e_d = [.54]; 
% e_i = [.91];
% e_c = [1.0156]; 

% %Fumeron,2024
% e_d = [0.610]; 
% e_i = [1.09];
% e_c = [0.947]; 


e0 =  [0.77,0.628,0.581,0.540]'; %Cuitiño,2016
% e0 =  [0.847,0.681,0.625,0.577]'; %Brown,2024
% e0 =  [0.89,0.738,0.688,0.644]'; %Fumeron,2024

% Parámetros generales
phi_c = 34.4; % Ángulo de fricción crítico (grados)
phi_c_rad = deg2rad(phi_c); % Convertir a radianes

% Constante a de la ecuación de Von Wolffersdorff
a = ( sqrt(3) * (3 - sin(phi_c_rad))) / (2 * sqrt(2) * sin(phi_c_rad));

% -------------------------------
% 2. Procesar cada ensayo
% -------------------------------
alpha_vals = zeros(1, num_ensayos); % Para almacenar los valores de alpha

for i = 1:num_ensayos
    % Extraer datos del ensayo actual
    idx = ensayo_indices(i, :);
    eps1 = data(:, idx(1)); % Deformación axial
    q = data(:, idx(2));    % Esfuerzo desviador
    p = data(:, idx(3));    % Presión media
    epsv = data(:, idx(4)); % Deformación volumétrica

    % Filtrar datos no nulos
    valid = ~isnan(eps1) & ~isnan(q) & ~isnan(p) & ~isnan(epsv);
    eps1 = eps1(valid);
    q = q(valid);
    p = p(valid);
    epsv = epsv(valid);

    % -------------------------------
    % Cálculo de parámetros intermedios
    % -------------------------------
    % Relación eta
    eta = q./ p;

    % Ángulo de fricción movil (phi_mob)
    phi_mob = asin(3 * eta ./ (6 + eta));

    % Ángulo de fricción máximo (phi_p)
    phi_p = max(phi_mob);
    phi_p_grados = rad2deg(phi_p); % Convertir a grados si es necesario

    % Relación de esfuerzos máximos (K_P)
    K_P = (1 + sin(phi_p))./ (1 - sin(phi_p));

    % Cálculo de A (ecuación 3.19)
    A = (a.^2 / (K_P + 2).^2) * (1 - (K_P * (4 - K_P) / (5 * K_P - 2)));

    % Cálculo de tan(nu_p) (ecuación 3.18)
    
%     tan_nup = ( (2*(K_P-4+5*A*K_P.^2 - 2*A*K_P))./ ...
%                  ((5*K_P-2)*(1+2*A)) ) -1;
    tan_nup = 2*((K_P - 4) + A * K_P * (5 * K_P - 2)) ./ ...
          ((5 * K_P - 2) * (1 + 2 * A)) - 1;
      
    r_e = (e0(i) - e_d) ./ (e_c - e_d); % Relación de densidad relativa

    % Calcular alpha (ecuación 3.14)
    alpha = (1 / log(r_e)) * log(...
            (6 * ((K_P + 2).^2 + a.^2 * K_P * (K_P - 1 - tan_nup))) / ...
            (a * (5 * K_P - 2) * (K_P + 2) * sqrt(4 + 2 * (1 + tan_nup).^2)));

    % Guardar valor de alpha
    alpha_vals(i) = alpha;
end

% -------------------------------
% 3. Visualizar resultados
% -------------------------------
fprintf('Valores de alpha para los cuatro ensayos:\n');
disp(alpha_vals);

% Graficar resultados
figure;
plot(alpha_vals,'o','linewidth',1.5);
xlabel('Ensayo')
ylabel('Valor de \alpha');
title('Calibración del parámetro \alpha para cada ensayo');
grid on
%% Script para calcular el parámetro alpha para múltiples ensayos triaxiales 
%s3 VARIATION
clear all;clc

% -------------------------------
% 1. Importar datos desde Excel
% -------------------------------
filename = 'data_alpha.xlsm'; % Nombre del archivo Excel
sheetName = 'Hoja2'; % Nombre de la hoja
data = readmatrix(filename, 'Sheet', sheetName);

% Índices de las columnas para cada ensayo
% TMD DR VARIATION
ensayo_indices = [1, 2, 3, 4;
                  6, 7, 8, 9;
                 11, 12, 13, 14; 
                 16, 17, 18, 19];
% % TMD S3 VARIATION
% ensayo_indices = [1, 2, 3, 4;
%                   6, 7, 8, 9;
%                  11, 12, 13, 14];

             
num_ensayos = size(ensayo_indices, 1);

% %Cuitiño,2016
% e_d = [0.508]; 
% e_i = [.947];
% e_c = [0.824]; 

% Rozas,2024
% e_d = [.54]; 
% e_i = [.91];
% e_c = [1.0156]; 
e_d = [.67]; 
e_i = [1.212];
e_c = [1.054]; 

% %Fumeron,2024
% e_d = [0.610]; 
% e_i = [1.09];
% e_c = [0.947]; 


% e0 =  [0.77,0.628,0.581,0.540]'; %Cuitiño,2016
e0 =  [0.847,0.681,0.625,0.577]'; %Brown,2024
% e0 =  [0.89,0.738,0.688,0.644]'; %Fumeron,2024

% Parámetros generales
phi_c = 34.4; % Ángulo de fricción crítico (grados)
phi_c_rad = deg2rad(phi_c); % Convertir a radianes

% Constante a de la ecuación de Von Wolffersdorff
a = ( sqrt(3) * (3 - sin(phi_c_rad))) / (2 * sqrt(2) * sin(phi_c_rad));

% -------------------------------
% 2. Procesar cada ensayo
% -------------------------------
alpha_vals = zeros(1, num_ensayos); % Para almacenar los valores de alpha

for i = 1:num_ensayos
    % Extraer datos del ensayo actual
    idx = ensayo_indices(i, :);
    eps1 = data(:, idx(1)); % Deformación axial
    q = data(:, idx(2));    % Esfuerzo desviador
    p = data(:, idx(3));    % Presión media
    epsv = data(:, idx(4)); % Deformación volumétrica

    % Filtrar datos no nulos
    valid = ~isnan(eps1) & ~isnan(q) & ~isnan(p) & ~isnan(epsv);
    eps1 = eps1(valid);
    q = q(valid);
    p = p(valid);
    epsv = epsv(valid);

    % -------------------------------
    % Cálculo de parámetros intermedios
    % -------------------------------
    % Relación eta
    eta = q./ p;

    % Ángulo de fricción movil (phi_mob)
    phi_mob = asin(3 * eta ./ (6 + eta));

    % Ángulo de fricción máximo (phi_p)
    phi_p = max(phi_mob);
    phi_p_grados = rad2deg(phi_p); % Convertir a grados si es necesario

    % Relación de esfuerzos máximos (K_P)
    K_P = (1 + sin(phi_p))./ (1 - sin(phi_p));

    % Cálculo de A (ecuación 3.19)
    A = (a.^2 / (K_P + 2).^2) * (1 - (K_P * (4 - K_P) / (5 * K_P - 2)));

    % Cálculo de tan(nu_p) (ecuación 3.18)
    
%     tan_nup = ( (2*(K_P-4+5*A*K_P.^2 - 2*A*K_P))./ ...
%                  ((5*K_P-2)*(1+2*A)) ) -1;
    tan_nup = 2*((K_P - 4) + A * K_P * (5 * K_P - 2)) ./ ...
          ((5 * K_P - 2) * (1 + 2 * A)) - 1;
      
    r_e = (e0 - e_d) ./ (e_c - e_d); % Relación de densidad relativa

    % Calcular alpha (ecuación 3.14)
    alpha = (1 / log(r_e)) * log(...
            (6 * ((K_P + 2).^2 + a.^2 * K_P * (K_P - 1 - tan_nup))) / ...
            (a * (5 * K_P - 2) * (K_P + 2) * sqrt(4 + 2 * (1 + tan_nup).^2)));

    % Guardar valor de alpha
    alpha_vals(i) = alpha;
end

% -------------------------------
% 3. Visualizar resultados
% -------------------------------
fprintf('Valores de alpha para los cuatro ensayos:\n');
disp(alpha_vals);

% Graficar resultados
figure;
plot(alpha_vals,'o','linewidth',1.5);
xlabel('Ensayo')
ylabel('Valor de \alpha');
title('Calibración del parámetro \alpha para cada ensayo');
grid on;
%%

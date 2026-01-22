clear all
clc

e23 =[0.8458	0.8425	0.8394	0.8374	0.8360	0.8326	0.8304	0.8270	0.8243	0.8222	0.8202	0.8181	0.8153	0.8133	0.8115	0.8092	0.8069	0.8032	0.8004	0.7960	0.7912	0.7870	0.7834	0.7800];
e26 =[0.841	0.838	0.834	0.831	0.830	0.827	0.825	0.821	0.818	0.817	0.814	0.812	0.810	0.808	0.806	0.804	0.801	0.798	0.795	0.791	0.787	0.783	0.779	0.776];

p =[1.02	1.41	2.18	2.96	3.74	5.69	7.63	11.53	16.48	21.44	28.95	36.46	50.50	64.54	78.58	98.09	122.45	161.45	200.46	262.07	340.09	418.07	496.11	574.13]';
e = [0.86584	0.86333	0.86093	0.85916	0.85801	0.85508	0.85341	0.85080	0.84819	0.84621	0.84381	0.84203	0.83890	0.83640	0.83431	0.83159	0.82888	0.82491	0.82136	0.81635	0.81113	0.80675	0.80226	0.79850]';
lnp = log(p);

Cc = diff(e)./diff(log(p));
% p = (1-(2/3)*sin(35))*pmean;
%%
% Datos experimentales
p = [1.02, 1.41, 2.18, 2.96, 3.74, 5.69, 7.63, 11.53, 16.48, 21.44, ...
     28.95, 36.46, 50.50, 64.54, 78.58, 98.09, 122.45, 161.45, 200.46, ...
     262.07, 340.09, 418.07, 496.11, 574.13]';
e = [0.86584, 0.86333, 0.86093, 0.85916, 0.85801, 0.85508, 0.85341, ...
     0.85080, 0.84819, 0.84621, 0.84381, 0.84203, 0.83890, 0.83640, ...
     0.83431, 0.83159, 0.82888, 0.82491, 0.82136, 0.81635, 0.81113, ...
     0.80675, 0.80226, 0.79850]';

% Definir el modelo con la exponencial
modelo = @(params, p) params(1) * exp(-(3 * p / params(2)).^params(3));

% Valores iniciales para los parámetros [e0, hs, n]
param_init = [0.69, 1000000, 0.23];  % Ajustar según el conocimiento inicial

% Ejecutar lsqcurvefit
[param_ajustados, resnorm] = lsqcurvefit(modelo, param_init, p, e);

% Mostrar resultados
disp('Parámetros ajustados:');
disp(['e0: ', num2str(param_ajustados(1))]);
disp(['hs: ', num2str(param_ajustados(2))]);
disp(['n: ', num2str(param_ajustados(3))]);

% Graficar el ajuste
figure;
semilogx(p, e, 'bo', 'MarkerFaceColor', 'b');  % Datos experimentales
hold on;
p_fit = linspace(min(p), max(p), 100);
e_fit = modelo(param_ajustados, p_fit);
semilogx(p_fit, e_fit, 'r-', 'LineWidth', 1.5);  % Curva ajustada
hold off;
xlabel('Presión media p');
ylabel('Relación de vacíos e');
title('Ajuste de la curva de compresión normal');
legend('Datos experimentales', 'Modelo ajustado');
grid on;

%%

clear all
clc

e = [0.91, 0.91, 0.90,0.89,0.88,0.87,0.85]';
p = [12.5,25,50,100,200,400,800]';

semilogx(p,e,'o')



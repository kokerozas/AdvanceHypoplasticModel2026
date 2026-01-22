clear all
clc
x = [linspace(0,1,100)];

f2 = 1 - exp(-5*x);  % saturación rápida hacia 1
f3 = 1 ./ (1 + exp(-10*(x - 0.5)));  % tipo sigmoide
f1 = x;   % crecimiento lineal

figure; hold on; grid on;
plot(x, f1, 'r', 'DisplayName', 'Lineal');
plot(x, f2, 'b', 'DisplayName', 'Exponencial saturada');
plot(x, f3, 'g', 'DisplayName', 'Sigmoide');
legend;
xlabel('x (variable normalizada)');
ylabel('f(x)');
title('Prototipos de ecuaciones entre 0 y 1');
%%
dx = x(2) - x(1);        % paso
df1 = diff(f1)/dx;
df2 = diff(f2)/dx;
df3 = diff(f3)/dx;

x_mid = x(1:end-1);      % vector para derivadas

figure; hold on; grid on;
plot(x_mid, df1, 'r', 'DisplayName', 'df1/dx');
plot(x_mid, df2, 'b', 'DisplayName', 'df2/dx');
plot(x_mid, df3, 'g', 'DisplayName', 'df3/dx');
legend;
xlabel('x');
ylabel('df/dx');
title('Derivadas de las funciones prototipo');
%%
chi0 = 2;
A = 1.5;

chi1 = chi0 + A * f1;
chi2 = chi0 + A * f2;
chi3 = chi0 + A * f3;

figure; hold on; grid on;
plot(x, chi1, 'r', 'DisplayName', 'Chi lineal');
plot(x, chi2, 'b', 'DisplayName', 'Chi exp');
plot(x, chi3, 'g', 'DisplayName', 'Chi sigmoide');
legend;
xlabel('x');
ylabel('\chi(x)');
title('Evolución de \chi con distintas funciones f(x)');

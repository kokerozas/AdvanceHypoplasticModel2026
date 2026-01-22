    %% Driver to test liquefaction model and modification of I_d density factor
clear all; clc;

% input_data_biobio;
% init_state_biobio;
% % 
input_data_karlsuhe;
init_state_karlsuhe;

% init_data = load('D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2025\OED - oedometric test - karlsuhe\OED_karlsuhe.txt');
% init_data = load('D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2023\TCUE - strain cycles large amplitudes\TCUE_karlsuhe.txt');
init_data = load('D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2023\TCU - cyclic triaxial test karlsuhe\TCUI_karlsuhe.txt');
% init_data = load('D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2025\TCU - cyclic triaxial test bio bio sand\TCU_biobiosand.txt');
% init_data = load('D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2023\TMU - triaxial-monotonic-undrained\TMU_fine_sand.txt');
% init_data = load('D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2023\TMU2 - triaxial monotonic undrained\TMU2_fine_sand.txt');
% init_data = load('D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2023\TMD - triaxial-monotonic-drained\TMD_fine_sand.txt');
% init_data = load('D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2025\TMD_biobio_cuitino/tmd_biobio_cuitino.txt');
% init_data = load('D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2025c\TCU - Hector Saldaña\TCUI_biobiosaldaña.txt');

kk = 1;
for i = 9%5:10:length(init_data)%9 %3:4:length(init_data) %5:10:length(init_data)  
%% Init state for all condition

%for OED compresion
% y0(7:9,1)      = init_data(i,3);
% y0(13,1)       = init_data(i,3);

%for constant amplitud strain 
% path_info(3,1) = init_data(i,5);
% y0(7:9,1)      = init_data(i,3);
% y0(13,1)       = init_data(i,1);

%for stress controlled test (TMD)
% y0(13,1)       = init_data(i,3);

%for cyclic condition
y0(7:9,1)      = init_data(i,2);
% y0(13,1)       = init_data(i,1);

global log_internal_vars log_index
log_index = 1;
log_internal_vars = struct( ...
    'fs', [], 'I_pe', [], 'chi', [], ...
    'd_strain', [], 'dZZ', [], ...
    'Se', [], 'FF', [], ...
    'Fd', [], 'dev_strain_rate_norm', [], ...
    'qdev', [], 'p', []  ,'pr', [], ...
    'l_dot', [] ...
);
%% Driver model 
tic
% [SS,EE,INV_S,INV_E,HARD] = updateModel(y0,parms,nspb,path_info);
[SS,EE,INV_S,INV_E,HARD] = updateModel_stress(y0,parms,nspb,path_info,25,25);

% xyplot_vmh; hold on; 

%% Import data
% for i = 1
% name = ['D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2025\OED - oedometric test - karlsuhe\OE' num2str(i) '.dat'];
% name = ['D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2023\TCUE - strain cycles large amplitudes\TCUE' num2str(i) '.dat'];
name = ['D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2023\TCU - cyclic triaxial test karlsuhe\TCUI' num2str(i) '.dat'];
% name = ['D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2025\TCU - cyclic triaxial test bio bio sand\cyc40_r' num2str(i) '.dat'];
% name = ['D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2023\TMU - triaxial-monotonic-undrained\TMU' num2str(i) '.dat'];
% name = ['D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2023\TMU2 - triaxial monotonic undrained\TMU' num2str(i) '.dat'];
% name = ['D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2023\TMD - triaxial-monotonic-drained\TMD' num2str(i) '.dat'];
% name = ['D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2025\TMD_biobio_cuitino\DR' num2str(i) '.dat'];
% name = ['D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2025c\TCU - Hector Saldaña\CSR' num2str(i) '.dat'];
a = importdata(name)
%%
% xyplot_vmh_tmd; hold on;
% xyplot_vmh_oed; hold on;

%oed 
% s1 = a.data(:,1);
% eps1 = a.data(:,2)/100;
% e = a.data(:,3);
% subplot(1,2,1)
% semilogx(s1,e,'o','color',[0.6 0.6 0.6]);hold on
% semilogx(s1,e,'bo');

% TMD (BIO-BIO SAND)
% eps1 = a.data(:,1)/100;
% epsv = a.data(:,4)/100;
% pexp = a.data(:,3);
% qexp = a.data(:,2);

% subplot(1,2,1)
% plot(eps1(1:end),qexp(1:end),'-','color',[0.00,0.45,0.74]);hold on
% % xlim([-0.15 105]); 
% subplot(1,2,2)
% plot(eps1,epsv,'-','color',[0.00,0.45,0.74]);hold on

%%  TCU - Rozas et al., 2024
% xyplot_vmh_cyc; hold on;
% eps1exp = a(:,2)/100;
% u = a(:,4)/100;
% pexp = a(:,3);
% qexp = a(:,4);

% %% FOR UNDRAINED CYCLIC TRIAXIAL TEST - KARLSUHE SAND PLOT
eps1 = a.data(:,9)/100;
p = a.data(:,7); 
q = a.data(:,8);
pp = a.data(:,4)/1;
xyplot_vmh_cyc; hold on;

subplot(2,2,1)
% plot(eps1,q,'k-');hold on
plot(p(1:end),q(1:end),'k-');hold on

subplot(2,2,2)
% plot(eps1,epsv,'k-');hold on
plot(eps1(1:end),q(1:end),'k-');hold on
xlim([-0.08, 0.08])

%$ TCU - rozas 
% xyplot_vmh_cyc; hold on;
% eps1exp = a(:,2)/100;
% u = a(:,4)/100;
% pexp = a(:,3);
% qexp = a(:,4);

%% TCU BIOBIOSAND HECTOR SALDAÑA
% eps1exp = a(:,1)/100;
% u = a(:,4);
% pexp = a(:,2);
% qexp = a(:,3);
% 
% figure(1); clf(1);
% subplot(2,2,1)
% plot(pexp,qexp,'color',[0.8 0.8 0.8]);hold on
% xline(0,'--','color',[0.7 0.7 0.7]); yline(0,'--','color',[0.7 0.7 0.7]);
% xlim([-10, 120]);
% yticks(-100:50:100); ylim([-100,100]);
% % yticks(-60:20:60); ylim([-60,60]);
% 
% xlabel('$p$ (kPa)', 'interpreter', 'latex');ylabel('$ q$ (kPa)','interpreter','latex')
% text(0.01, 0.98,...
%     {'{Experiment data}', ...
%      ''}, ...
%     'Units', 'normalized', ...
%     'Interpreter', 'latex', ...
%     'FontSize', 10, ...
%     'HorizontalAlignment', 'left', ...
%     'VerticalAlignment', 'top')
% 
% text(-0.15, 1.05, '(A)', ...
%     'Units', 'normalized', ...
%     'FontWeight', 'bold', ...
%     'FontSize', 10, ...
%     'Interpreter', 'Latex', ...
%     'HorizontalAlignment', 'left', ...
%     'VerticalAlignment', 'top', ...
%     'Clipping', 'off')
% 
% subplot(2,2,2)
% plot(eps1exp,qexp,'color',[0.8 0.8 0.8]);hold on
% xline(0,'--','color',[0.7 0.7 0.7]); yline(0,'--','color',[0.7 0.7 0.7]);
% xlabel('$\varepsilon_1$ \, (-)', 'Interpreter', 'latex'); ylabel('$ q$ (kPa)','interpreter','latex')
% xlim([-0.06, 0.06]);
% % yticks(-60:20:60); ylim([-60,60]);ylabel('$q$ (kPa)','interpreter','latex')
% yticks(-100:50:100); ylim([-100,100]);
% text(0.01, 0.98,...
%     {'{Experiment data}', ...
%      ''}, ...
%     'Units', 'normalized', ...
%     'Interpreter', 'latex', ...
%     'FontSize', 10, ...
%     'HorizontalAlignment', 'left', ...
%     'VerticalAlignment', 'top')
% 
% text(-0.15, 1.05, '(B)', ...
%     'Units', 'normalized', ...
%     'FontWeight', 'bold', ...
%     'FontSize', 10, ...
%     'Interpreter', 'Latex', ...
%     'HorizontalAlignment', 'left', ...
%     'VerticalAlignment', 'top', ...
%     'Clipping', 'off')
% 
% xyplot_vmh_cyc; hold on;
% 
% set(gcf,'Name','2025c','NumberTitle','off');
% hold on;
% 
% 
% set(gcf, 'Units', 'centimeters')
% set(gcf, 'Position', [0 0 26 16])  % Posición y tamaño de ventana
%% EXPORT TMD
% [~, pos] = unique(eps1);
% epsv_in(:,kk) = interp1(eps1(pos), epsv(pos), EE(:,3));
% q_in(:,kk)    = interp1(eps1(pos), qexp(pos)   , EE(:,3));
% 
% writematrix(epsv_in,'epsv_in_biobio.txt');
% writematrix(q_in,'q_in_biobio.txt');
%% Export cyclic data 
% Inputs: vector_model y nuevo tamaño deseado
% n_model = 54418; % data points to output interpolate
% % === Crear nuevo vector base de interpolación ===
% x_old = linspace(0, 1, length(q));
% x_new = linspace(0, 1, n_model);
% % Export vector data
% eps1_in_cyc8 = interp1(x_old, eps1, x_new, 'spline');
% p_in_cyc8 = interp1(x_old, p, x_new, 'spline');
% q_in_cyc8 = interp1(x_old, q, x_new, 'spline');
% % Save data
% writematrix(eps1_in_cyc8,'eps1_in_cyc8.txt');
% writematrix(p_in_cyc8,'p_in_cyc8.txt');
% writematrix(q_in_cyc8,'q_in_cyc8.txt');
%% pause()
kk = 1 + kk;
tiempo_compute = toc;
fprintf('El tiempo de computo es (s): %.4f\n', tiempo_compute);
fprintf('El tiempo de computo es (min): %.4f\n', tiempo_compute/60);
end
%%
%plot(log_internal_vars.chi)
%% Calibration of \eta_2 and C_{\chi} in pre liquefaction stage
% ---CICLOS EXPERIMENTALES PRE LICUACIÓN---
nq = length(qexp);
ciclovector = zeros(nq,1);

% Detectar cambios de signo (pquntos de cruce por cero)

sign_q = diff(sign(qexp));
crossidx = find(sign_q ~= 0);

crossidx = [1; crossidx(:); nq];  % incluir extremos

% Contar y asignar ciclos
cicloval = 0;
for i = 1:length(crossidx)-1
    idxrange = crossidx(i):crossidx(i+1);
    ciclovector(idxrange) = linspace(cicloval, cicloval+0.5, length(idxrange));
    cicloval = cicloval + 0.5;
    
    
end
nciclosestimados = cicloval;

fprintf('Número de ciclos experimentales pre licuación: %.0f\n', nciclosestimados);

% --- PRESIÓN DE POROS ACUMULADA ---
p0 = 100; % confinamiento inicial [kPa]
pexp1 = pexp(1:end);     % presión total experimental (pre licuación)
uexp = p0 - pexp1;
% pmodel = log_internal_vars.p(:);

pmodel = (SS(:,3)+2*SS(:,1))/3%INV_S(:,1);
umodel = p0 - pmodel;

% Suavizado
umodel_smooth = movmean(umodel,100);
uexp_smooth   = movmean(uexp, 1);

%--- CICLOS DEL MODELO (sin cambios) ---
% strain = log_internal_vars.eps_a(:);
% qdev = log_internal_vars.qdev(:);
% strain = eps1(:);
strain = EE(:,3);
% qdev = SS(:,3)-SS(:,1);

n_points = length(strain);
ciclo_vector = zeros(n_points, 1);
sign_model = diff(sign(strain));
cross_idx = find(sign_model ~= 0);
cross_idx = [1; cross_idx(:); n_points];
ciclo_val = 0;
for i = 1:length(cross_idx)-1
    idx_range = cross_idx(i):cross_idx(i+1);
    ciclo_vector(idx_range) = linspace(ciclo_val, ciclo_val+0.5, length(idx_range));
    ciclo_val = ciclo_val + 0.5;
end
n_ciclos_estimados = ciclo_val;
fprintf('Número de ciclos modelo total: %.0f\n', n_ciclos_estimados);

% --- GRAFICAR p_w vs N ---
% subplot(3,4,12)
subplot(1,1,1)
plot(ciclo_vector, umodel_smooth, '--','linewidth',0.5); hold on
plot(ciclovector, uexp_smooth,'color',[0.75 0.75 0.75]); hold on
 
% xticks(0:20:60); xlim([0,60]);%tcui9
yticks(0:20:100); ylim([0,100]);
xlabel('$N$ (-)','interpreter','latex')
ylabel('$p_w^{acc}$ (kPa)','interpreter','latex')

% set(gcf, 'Units', 'centimeters')
% set(gcf, 'Position', [5 2 25 15])  % Posición y tamaño de ventana
%%
% === Cargar parámetros e inicializar ===
input_data_biobio;
init_state_biobio;

% === Combinaciones de parámetros ===
mr  = 2.2;
chi = 5.0;
br  = 0.15;
A   = combvec(mr, chi, br)';   % [mr, chi, betar]

% === Nivel de confinamiento (kPa) ===
sigv = [25 50 100 200 400];
phi  = 33.1;                    % Ángulo de fricción (°)
k0   = 1 - sind(phi);           % Coef. de presiones en reposo

% === Inicializar resultados ===
g0 = zeros(length(sigv), size(A,1));

% === Loop principal ===
for j = 1:length(sigv)
    y0(7:9,1) = [sigv(j); k0*sigv(j); k0*sigv(j)];
    
    for ik = 1:size(A,1)
        % === Parámetros actualizados ===
%         parms(12)  = 2.2;                    % mR
%         parms(13)  = 2.2 * A(ik,1);          % Rparam
%         parms(15)  = A(ik,3);                % betar
%         parms(16)  = A(ik,2);                % chi

        % === Simulación ===
%         [SS, EE, ~, ~, ~] = updateModel_stress(y0, parms, nspb, path_info, 25, 25);
        [SS, EE, ~, ~, ~] = updateModel(y0, parms, nspb, path_info);
        q     = SS(:,3) - SS(:,1);      % q = σ11 - σ33
        eps1  = EE(:,3);                % ε1 axial

        % === Cálculo rigidez tangente secante ===
        G = NaN(length(q)-1, 1);
        for i = 1:length(q)-1
            de = eps1(i+1) - eps1(i);
            dq = q(i+1) - q(i);
            if abs(de) > 1e-6
                G(i) = dq / de;
            end
        end

        g0(j,ik) = max(G);  % guardar Gmáx para esta combinación

        % === Plot rigidez ===
        figure(3)
        semilogx(abs(eps1(2:end)), G, 'DisplayName', ...
            ['\sigma_v = ' num2str(sigv(j)) ', m_r=' num2str(A(ik,1)) ', \chi=' num2str(A(ik,2))]); 
        hold on; grid on;
        xlabel('\epsilon_1'); ylabel('G = \Delta q / \Delta \epsilon_1');
        title('Rigidez tangente axial G vs \epsilon_1');
    end
end

legend('Location', 'best');
%%
subplot(2,4,5)
plot(log_internal_vars.Fd,'r-'); hold on
subplot(2,4,6)
plot(log_internal_vars.ag,'r-'); hold on
subplot(2,4,7)
plot(log_internal_vars.a01,'r-'); hold on
subplot(2,4,8)
plot(log_internal_vars.chi,'b-'); hold on


%% II PART (DENSE MATERIALS)
% --- CICLOS EXPERIMENTALES PRE LICUACIÓN ---
exp = qexp(1:end);  % señal de esfuerzo desviado pre licuación
nq = length(qexp);
ciclovector = zeros(nq,1);


% Detectar cambios de signo (puntos de cruce por cero)
sign_q = diff(sign(qexp));
crossidx = find(sign_q ~= 0);
crossidx = [1; crossidx(:); nq];  % incluir extremos

% Contar y asignar ciclos
cicloval = 0;
for i = 1:length(crossidx)-1
    idxrange = crossidx(i):crossidx(i+1);
    ciclovector(idxrange) = linspace(cicloval, cicloval+0.5, length(idxrange));
    cicloval = cicloval + 0.5;
    
end
nciclosestimados = cicloval;

fprintf('Número de ciclos experimentales pre licuación: %.0f\n', nciclosestimados);

% --- CORTAR PRIMEROS 10 CICLOS EXPERIMENTALES ---
Ndrop = 19;                                           % ciclos a eliminar
if max(ciclovector) < Ndrop
    warning('Solo hay %.1f ciclos; no se recorta.', max(ciclovector));
else
    mask = ciclovector >= Ndrop;                      % conservar desde N>=10
    % re-referenciar N para que el nuevo inicio sea N=0
    ciclovector = ciclovector(mask) - Ndrop;

    % recortar series EXPERIMENTALES asociadas
    qexp = qexp(mask);
    pexp = pexp(mask);
    exp  = exp(mask);   % si la usas
end

% --- PRESIÓN DE POROS ACUMULADA ---
p0 = 100; % confinamiento inicial [kPa]
pexp1 = pexp(1:end);     % presión total experimental (pre licuación)
uexp = p0 - pexp1;
% pmodel = log_internal_vars.p(:);

pmodel = (SS(:,3)+2*SS(:,1))/3%INV_S(:,1);
umodel = p0 - pmodel;

% Suavizado
umodel_smooth = movmean(umodel,100);
uexp_smooth   = movmean(uexp, 1);

%--- CICLOS DEL MODELO (sin cambios) ---
% strain = log_internal_vars.eps_a(:);
% qdev = log_internal_vars.qdev(:);
% strain = eps1(:);
strain = EE(:,3);
% qdev = SS(:,3)-SS(:,1);

n_points = length(strain);
ciclo_vector = zeros(n_points, 1);
sign_model = diff(sign(strain));
cross_idx = find(sign_model ~= 0);
cross_idx = [1; cross_idx(:); n_points];
ciclo_val = 0;
for i = 1:length(cross_idx)-1
    idx_range = cross_idx(i):cross_idx(i+1);
    ciclo_vector(idx_range) = linspace(ciclo_val, ciclo_val+0.5, length(idx_range));
    ciclo_val = ciclo_val + 0.5;
end
n_ciclos_estimados = ciclo_val;
fprintf('Número de ciclos modelo total: %.0f\n', n_ciclos_estimados);

% --- GRAFICAR p_w vs N ---
subplot(3,4,12)
plot(ciclo_vector, umodel_smooth, 'r-','linewidth',0.5); hold on
plot(ciclovector, uexp_smooth, 'k-','color',[0.75 0.75 0.75]); hold on
% xline(10);hold on
% xticks(0:10:60); xlim([0,60]); %tcui 18
xticks(0:5:15); xlim([0,15]); %tcui 18
% xticks(0:20:100); xlim([0,100]);
yticks(0:10:70); ylim([0,70]);
xlabel('$N$ (-)','interpreter','latex')
ylabel('$p_w^{acc}$ (kPa)','interpreter','latex')
legend('Modelo','Experimental')
set(gca, 'LineWidth', 0.7)
% text(-0.18, 1.05, '(L)', ...
%     'Units', 'normalized', ...
%     'FontWeight', 'bold', ...
%     'FontSize', 12, ...
%     'Interpreter', 'Latex', ...
%     'HorizontalAlignment', 'left', ...
%     'VerticalAlignment', 'top', ...
%     'Clipping', 'off')

% set(gcf, 'Units', 'centimeters')
% set(gcf, 'Position', [5 2 25 15])  % Posición y tamaño de ventana
%% === Inicialización de la calibración: TCUI9 (eps1–q) ===
clear; close all; clc
clear functions; rehash

% 1) Cargar datos experimentales desde TXT (columnas)
p_exp    = readmatrix('p_in_cyc9.txt');    p_exp    = p_exp(:);
q_exp    = readmatrix('q_in_cyc9.txt');    q_exp    = q_exp(:);
eps1_exp = readmatrix('eps1_in_cyc9.txt'); eps1_exp = eps1_exp(:);

% 2) Preparar modelo base (NO cambiar driver)
input_data_karlsuhe;    % crea: parms, nspb, path_info   (¡sin "b"!)
init_state_karlsuhe;    % crea: y0                       (¡sin "b"!)


% 3) Semilla y función objetivo (coherentes con la fitness)
x0 = [1000, 10.45, 0.2];   % [cz, zmax, lambda1]
fun = @(x) fitnessfunction_karlsuhe_cyc2(x, y0, parms, nspb, path_info, p_exp, q_exp, eps1_exp);


% 4A) fminsearch (rápido, sin restricciones; la fitness penaliza)
% opts = optimset('Display','iter', 'TolX',1e-3, 'TolFun',1e-4, ...
%                 'MaxFunEvals',10000, 'MaxIter',2000);
% [xopt,Jmin,exitflag,output] = fminsearch(fun, x0, opts);


% % 4B) (alternativa) fmincon con límites
lb = [500, 1.0, 0.10];  ub = [2000, 15.0, 1.0];
opts = optimoptions('fmincon','Display','iter','Algorithm','sqp', ...
        'MaxFunctionEvaluations',3e4,'MaxIterations',2000);
[xopt,Jmin,exitflag,output] = fmincon(fun, x0, [],[],[],[], lb, ub, [], opts);

% 5) Reporte
fprintf('\nÓPTIMO: cz=%.0f, zmax=%.2f, lambda1=%.3f | J=%.3e\n', xopt(1), xopt(2), xopt(3), Jmin);
disp(output)

%%
%% Improvment parameters
% alpha(ee)  = -0.482*x^18.63+0.6103;
% beta(ee)   = 1.25*x^16.63+6.344;
% rparam(ee) = (1.531*x^2 + -2.366*x + 0.9191) / (x^2 + -1.565*x + 0.6181);
% q_br(ee)   = 0.609*exp(-0.8411*x) + 1.278e-16*exp(35.28*x)
% En(ee)     =  0.1864 + 0.01375*cos(x*19.63) + -0.1*sin(x*19.63)

% x0 = [1.1757    5.1615    1.2428   19.3490    0.0979   16.8263    0.8719]
x1 = [0.5  0   0.5  0.1   0     10  ];
x2 = [1.2  12  1.7  3.0   0.15  35  ];
% A  = [0 0 0 0 0 0  1 -1  0;
%       0 0 0 0 0 0 -1  0  1];
% b  = [0;0];
%% Monotonic optimization block
% input_data_biobio;
% init_state_biobio;

input_data_karlsuhe;
init_state_karlsuhe;

TAU    = load('q_in_biobio.txt');
DEV    = load('epsv_in_biobio.txt');
% TAU    = load('q_in.txt');
% DEV    = load('epsv_in.txt');

ObjFun = @(x) fitnessfunction_biobio(x,TAU,DEV,init_data);
% ObjFun = @(x) fitnessfunction_karlsuhe(x,TAU,DEV,init_data);
% options = optimset('TolFun',0.8,'TolX',1,'OutputFcn', @myoutput); %opt1
options = optimset('OutputFcn', @myoutput); %opt2
options.MaxFunEvals = 500;
options.MaxIter     = 500;
options.UseParallel = 1;


%    alpha   r_parm  beta    epsmon
% x0 = [0.77   1.7     0.12    14.5]
     
%    alpha   r_parm  beta    epsmon
x0 = [0.77   1.7    0.12    14.5]

[x , fval] = fminsearch(ObjFun,x0,options);
% [x , fval] = patternsearch(ObjFun,x0,[],[],[],[],x1,x2,[],options);
%% Shear moduls analysis

% === Cargar parámetros e inicializar ===
input_data_karlsuhe;
init_state_karlsuhe;

% === Combinaciones de parámetros ===
mr  = 2.2;
chi = 5.0;
br  = 0.15;
A   = combvec(mr, chi, br)';   % [mr, chi, betar]

% === Nivel de confinamiento (kPa) ===
sigv = [25 50 100 200 400];
phi  = 33.1;                    % Ángulo de fricción (°)
k0   = 1 - sind(phi);           % Coef. de presiones en reposo

% === Inicializar resultados ===
g0 = zeros(length(sigv), size(A,1));

% === Loop principal ===
for j = 1:length(sigv)
    y0(7:9,1) = [sigv(j); k0*sigv(j); k0*sigv(j)];
    
    for ik = 1:size(A,1)
        % === Parámetros actualizados ===
%         parms(12)  = 2.2;                    % mR
%         parms(13)  = 2.2 * A(ik,1);          % Rparam
%         parms(15)  = A(ik,3);                % betar
%         parms(16)  = A(ik,2);                % chi

        % === Simulación ===
        [SS, EE, ~, ~, ~] = updateModel_stress(y0, parms, nspb, path_info, 25, 25);
%         [SS, EE, ~, ~, ~] = updateModel(y0, parms, nspb, path_info);
        q     = SS(:,3) - SS(:,1);      % q = σ11 - σ33
        eps1  = EE(:,3);                % ε1 axial

        % === Cálculo rigidez tangente secante ===
        G = NaN(length(q)-1, 1);
        for i = 1:length(q)-1
            de = eps1(i+1) - eps1(i);
            dq = q(i+1) - q(i);
            if abs(de) > 1e-6
                G(i) = dq / de;
            end
        end

        g0(j,ik) = max(G);  % guardar Gmáx para esta combinación

        % === Plot rigidez ===
        figure(3)
        semilogx(abs(eps1(2:end)), G, 'DisplayName', ...
            ['\sigma_v = ' num2str(sigv(j)) ', m_r=' num2str(A(ik,1)) ', \chi=' num2str(A(ik,2))]); 
        hold on; grid on;
        xlabel('\epsilon_1'); ylabel('G = \Delta q / \Delta \epsilon_1');
        title('Rigidez tangente axial G vs \epsilon_1');
    end
end

legend('Location', 'best');



%%

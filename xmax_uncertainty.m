%% ===== chimax_uncertainty.m =====
function xmax_uncertainty
% Ley: log(chi) = β0 + β1 log(D_r) + β2 log(CSR)
% (A) chi vs D_r: franjas-ley por CSR + IC95 (todas las pruebas) + medianas (LW=0.1)
% (B) chi vs CSR: franjas-ley por D_r + IC95 (todas las pruebas, con clipping) + medianas (LW=0.1)
% (C) leave-one-out con banda operativa ±5 %

clear; close all; clc; rng(1)


%% -------------------- Estilo / opciones --------------------
tol        = 0.05;                       % banda operativa (panel C)
fontName   = 'Times New Roman'; fontSize = 11;
set(0,'DefaultAxesFontName',fontName,'DefaultAxesFontSize',fontSize);
set(0,'DefaultAxesTickDir','out','DefaultAxesLineWidth',0.5);

icAlpha    = 0.20;                       % transparencia IC95
rangeAlpha = 0.14;                       % transparencia franjas-ley
colMustard = [0.89 0.67 0.16];           % ley CSR o D_r (mostaza)
colGreen   = [0.00 0.62 0.45];           % ley CSR o D_r (verde/teal)

CSR_bands  = [0.20 0.30; 0.30 0.50];     % franjas-ley en Panel A
Dr_bands   = [0.60 0.70; 0.70 0.80];     % franjas-ley en Panel B

%% -------------------- Datos (χ_max calibrado) --------------------
tests   = ["tcui8","tcui9","tcui10","tcui20","tcui19","tcui18"]';
Dr_pts  = [0.62  0.68  0.61  0.79  0.78  0.78]';
CSR_pts = [0.20  0.25  0.30  0.50  0.40  0.30]';
chi_cal = [13.7 16.7 26.8 33.53 26.8 16.0]';           % χ_max medido/calibrado

nums    = regexp(tests,'\d+','match','once');
labels  = "TCUI " + string(nums);

%% -------------------- Diseño en log y OLS --------------------
Y = log(chi_cal);
X = [ones(size(Dr_pts))  log(Dr_pts)  log(CSR_pts)];
[n,p] = size(X);                          % n=6, p=3

bOLS  = (X'*X)\(X'*Y);                    % [β0; β1; β2]
yhat  = X*bOLS;
res   = Y - yhat;
H     = X/(X'*X)*X';  h = diag(H);

a = exp(bOLS(1)); b = bOLS(2); c = bOLS(3);
chifun = @(Dr,CSR) a.*(Dr.^b).*(CSR.^c);

% LOO exacto (panel C)
yhat_LOO = Y - res./(1 - h);
pred_LOO = exp(yhat_LOO);
rel      = (pred_LOO - chi_cal)./chi_cal;
MAPE     = 100*mean(abs(rel));
RMSErel  = 100*sqrt(mean(rel.^2));
biasRel  = 100*mean(rel);
%% ========= SENSIBILIDAD A chi_cal(TCUI20) — VISUALIZACIÓN =========
% % Este bloque NO altera los gráficos previos; genera una figura nueva.
% 
% % --- 1) Detectar el vector observado (chi_cal o z_cal) y la etiqueta ---
% if exist('chi_cal','var')
%     yobs_name = '$\chi_{\max}$';
%     yobs_base = chi_cal(:);
% elseif exist('chi_exp','var')
%     yobs_name = '$\chi_{\max}$';
%     yobs_base = chi_exp(:);
% else
%     yobs_name = '$z_{\max}$';
%     yobs_base = z_cal(:);
% end
% n = numel(yobs_base);
% 
% % --- 2) Índice de TCUI20 y valor modificado ---
% ix20      = find(strcmpi(tests,'tcui20'), 1);
% if isempty(ix20)
%     error('No se encontró TCUI20 en "tests".');
% end
% y_old_20  = yobs_base(ix20);   % valor original (p.ej. 31.5)
% y_new_20  = 35;              % <-- cambia aquí tu escenario
% 
% % --- 3) Baseline: re-calcular LOO y métricas (por si no están en workspace) ---
% %     Usamos el mismo X de tu ajuste base y recomputamos por claridad.
% Y_base     = log(yobs_base);
% b_base     = (X'*X)\(X'*Y_base);
% res_base   = Y_base - X*b_base;
% H          = X/(X'*X)*X';  h = diag(H);
% yhatLOO_b  = Y_base - res_base./(1 - h);         % LOO exacto (en log)
% predLOO_b  = exp(yhatLOO_b);                     % escala real
% er_b       = 100*(predLOO_b - yobs_base)./yobs_base;  % error relativo %
% 
% % --- 4) Nuevo: solo TCUI20 cambia (31.5 → 32.5) y re-estima ---
% Y_new            = Y_base;
% Y_new(ix20)      = log(y_new_20);
% b_new            = (X'*X)\(X'*Y_new);
% res_new          = Y_new - X*b_new;
% yhatLOO_new      = Y_new - res_new./(1 - h);
% predLOO_new      = exp(yhatLOO_new);
% er_new           = 100*(predLOO_new - yobs_base)./yobs_base;  % *mismo denominador* para comparar
% d_er             = er_new - er_b;                              % cambio por ensayo (puntos porcentuales)
% 
% % --- 5) Coeficientes (a,b,c) base vs nuevo ---
% a_b = exp(b_base(1));  b_b = b_base(2);  c_b = b_base(3);
% a_n = exp(b_new(1));   b_n = b_new(2);   c_n = b_new(3);
% 
% % --- 6) Gráfico: 1x2, errores y deltas ---
% figS = figure('Color','w','Units','centimeters','Position',[2 2 26 9]);
% tloS = tiledlayout(figS,1,2,'TileSpacing','compact','Padding','compact');
% 
% % ===== (S1) Errores relativos por ensayo (base vs nuevo) =====
% axS1 = nexttile(tloS,1); hold(axS1,'on'); box(axS1,'off')
% 
% % Banda ±tol (en %)
% xpos = 1:n;
% pBand = patch(axS1,[0.5 n+0.5 n+0.5 0.5], [-100*tol -100*tol 100*tol 100*tol], ...
%               [0.93 0.95 0.98], 'EdgeColor','none', 'FaceAlpha',0.65, ...
%               'DisplayName', sprintf('$\\pm %g\\,\\%%$', 100*tol));
% 
% % Segmentos base→nuevo y puntos
% for i=1:n
%     plot(axS1, [xpos(i) xpos(i)], [er_b(i) er_new(i)], ':', 'Color',[0.4 0.4 0.4], 'LineWidth',0.9, ...
%          'HandleVisibility','off');
% end
% pBase = scatter(axS1, xpos, er_b, 36, 'k','filled', 'DisplayName','Base');
% pNew  = scatter(axS1, xpos, er_new, 28, [0.00 0.45 0.74], 'filled', 'DisplayName','Nuevo');
% 
% % Resaltar TCUI20
% plot(axS1, [xpos(ix20) xpos(ix20)], [er_b(ix20) er_new(ix20)], '-', ...
%      'Color',[0.00 0.45 0.74], 'LineWidth',1.2, 'HandleVisibility','off');
% scatter(axS1, xpos(ix20), er_new(ix20), 48, [0.00 0.45 0.74], 'filled', 'MarkerEdgeColor','k', ...
%         'DisplayName','TCUI20');
% 
% % Ejes/rotulación
% xlim(axS1,[0.5 n+0.5]);
% yl1 = max(abs([er_b; er_new; 100*tol]))*1.15;  ylim(axS1,[-yl1 yl1])
% set(axS1,'XTick',xpos,'XTickLabel',upper(string(tests)));
% xtickangle(axS1,25)
% ylabel(axS1,'Relative error $e_r$ (\%)','Interpreter','latex')
% title(axS1, sprintf('Sensitivity at TCUI20: %s %g \\(\\to\\) %g', yobs_name, y_old_20, y_new_20), ...
%       'Interpreter','latex','FontSize',10)
% legend(axS1,[pBand pBase pNew], {'$\pm$ tolerance','Base','New'}, ...
%        'Interpreter','latex','Location','southoutside','Orientation','horizontal','Box','off');
% 
% % Métricas resumen (opcionales) en texto
% MAPE_b  = mean(abs(er_b));                 RMSE_b  = sqrt(mean(er_b.^2));                 Bias_b  = mean(er_b);
% MAPE_n  = mean(abs(er_new));               RMSE_n  = sqrt(mean(er_new.^2));               Bias_n  = mean(er_new);
% txt = sprintf(['Base: MAPE=%.2f%%, RMSE_{rel}=%.2f%%, Bias=%.2f%%\\n',...
%                'New : MAPE=%.2f%%, RMSE_{rel}=%.2f%%, Bias=%.2f%%'], ...
%                MAPE_b, RMSE_b, Bias_b, MAPE_n, RMSE_n, Bias_n);
% text(axS1, 0.55, -yl1+0.06*(2*yl1), txt, 'Interpreter','latex', 'FontSize',9, 'Color',[0.15 0.15 0.15]);
% 
% % ===== (S2) Cambio de error por ensayo (Δe_r) =====
% axS2 = nexttile(tloS,2); hold(axS2,'on'); box(axS2,'off')
% bBar = bar(axS2, xpos, d_er, 0.70, 'FaceColor',[0.70 0.70 0.70], 'EdgeColor','none', ...
%            'DisplayName','$\Delta e_r$');
% % resaltar TCUI20
% bBar.CData(ix20,:) = repmat([0.00 0.45 0.74],1,1);
% 
% yline(axS2,0,'k-','LineWidth',0.8,'HandleVisibility','off')
% xlim(axS2,[0.5 n+0.5]);
% yl2 = max(abs(d_er))*1.15; ylim(axS2,[-yl2 yl2])
% set(axS2,'XTick',xpos,'XTickLabel',upper(string(tests)));
% xtickangle(axS2,25)
% ylabel(axS2,'$\Delta e_r$ (p.p.)','Interpreter','latex')
% title(axS2,'Change in relative error (New − Base)','Interpreter','latex','FontSize',10)
% legend(axS2,bBar, {'$\Delta e_r$'}, 'Interpreter','latex','Location','southoutside', ...
%        'Orientation','horizontal','Box','off');
% 
% % --- 7) (Opcional) imprimir a consola el cambio en (a,b,c) ---
% fprintf('\nSensibilidad TCUI20: %s %.3f -> %.3f\n', erase(yobs_name,'$'), y_old_20, y_new_20);
% fprintf('Coeficientes base : a=%.3f, b=%.3f, c=%.3f\n', a_b, b_b, c_b);
% fprintf('Coeficientes nuevos: a=%.3f, b=%.3f, c=%.3f\n', a_n, b_n, c_n);

%% -------------------- Bootstrap de parámetros (wild Rademacher) --------------------
B = 5000;                                  % 2000–5000 ok para n=6
u_tilde = res ./ (1 - h);                  % residuo corregido
bBoot   = zeros(B,3);
XtX_inv = inv(X'*X);                       % eficiencia
for k = 1:B
    v    = sign(rand(n,1)-0.5);            % ±1
    Yb   = X*bOLS + u_tilde .* v;          % bootstrap en log
    bBoot(k,:) = XtX_inv*(X'*Yb);          % re-ajuste OLS
end

%% -------------------- Rangos y límites --------------------
Drv  = linspace(0.55,0.85,400);
CSRv = linspace(0.20,0.50,400);
[DRg, CSg] = ndgrid(Drv, CSRv);
Chi_all = chifun(DRg, CSg);
yl = [min([Chi_all(:); chi_cal])*0.95, max([Chi_all(:); chi_cal])*1.05];

%% -------------------- FIGURA 1×3 --------------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 28 10]);
tlo = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');

%% ================== (A) chi vs D_r ==================
ax1 = nexttile(tlo,1); hold(ax1,'on'); box(ax1,'off')

% Franjas-ley por CSR (al fondo)
y1 = chifun(Drv, CSR_bands(1,1)); y2 = chifun(Drv, CSR_bands(1,2));
hCSR1 = fill_between(Drv, min(y1,y2), max(y1,y2), colMustard, rangeAlpha, '$\mathrm{CSR}$ 0.20–0.30 (ley)');
y1 = chifun(Drv, CSR_bands(2,1)); y2 = chifun(Drv, CSR_bands(2,2));
hCSR2 = fill_between(Drv, min(y1,y2), max(y1,y2), colGreen,   rangeAlpha, '$\mathrm{CSR}$ 0.30–0.50 (ley)');
uistack([hCSR1 hCSR2],'bottom');

% Paleta en escala de grises para todas las IC95
grayLevels = linspace(0.85, 0.30, numel(tests));
cols = repmat(grayLevels', 1, 3);

% Bandas IC95 (todas las pruebas) + medianas (LW=0.1)
hIC_A = gobjects(1,1);
for j = 1:numel(tests)
    [Clo,~,Chi] = ci_curve_bootstrap(Drv, CSR_pts(j), bBoot, 'Dr', bOLS);
    h = fill_between(Drv, Clo, Chi, cols(j,:), icAlpha, 'IC95 (todas las pruebas)');
    if j==1, hIC_A = h; else, set(h,'HandleVisibility','off'); end

    % mediana IC95 (percentil 50) — grosor 0.1
    XcA   = [ones(numel(Drv),1), log(Drv(:)), log(CSR_pts(j))*ones(numel(Drv),1)];
    chiMedA = exp(prctile(XcA*bBoot.', 50, 2));
    plot(ax1, Drv, chiMedA, 'k-', 'LineWidth', 0.1, 'HandleVisibility','off');
end

% Datos
pDataA = scatter(ax1, Dr_pts, chi_cal, 44, 'k','filled','DisplayName','Datos calibrados');
for i=1:numel(chi_cal)
    text(ax1, Dr_pts(i)+0.006*range(Drv), chi_cal(i)+0.015*range(yl), labels(i), ...
         'FontName',fontName,'FontSize',7.5,'Clipping','on');
end

xlim(ax1,[min(Drv) max(Drv)]); ylim(ax1, yl)
xlabel(ax1,'$D_r$ (-)','Interpreter','latex'); 
ylabel(ax1,'$\chi_{\max}$ (-)','Interpreter','latex')
legend(ax1,[hCSR1 hCSR2 hIC_A pDataA], ...
       {'$\mathrm{CSR}$ 0.20–0.30 (ley)','$\mathrm{CSR}$ 0.30–0.50 (ley)','IC95 (todas)','Datos'}, ...
       'Interpreter','latex','Location','northwest','Box','off');

%% ================== (B) chi vs CSR (con clipping por D_r) ==================
ax2 = nexttile(tlo,2); hold(ax2,'on'); box(ax2,'off')

% Franjas-ley por D_r (al fondo)
y1 = chifun(Dr_bands(1,1)*ones(size(CSRv)), CSRv);
y2 = chifun(Dr_bands(1,2)*ones(size(CSRv)), CSRv);
hDR1 = fill_between(CSRv, min(y1,y2), max(y1,y2), colMustard, rangeAlpha, '$D_r$ 0.60–0.70 (ley)');
y1 = chifun(Dr_bands(2,1)*ones(size(CSRv)), CSRv);
y2 = chifun(Dr_bands(2,2)*ones(size(CSRv)), CSRv);
hDR2 = fill_between(CSRv, min(y1,y2), max(y1,y2), colGreen,   rangeAlpha, '$D_r$ 0.70–0.80 (ley)');
uistack([hDR1 hDR2],'bottom');

% Límites columna para clipping
DR1_lo = min( chifun(Dr_bands(1,1)*ones(size(CSRv)),CSRv), chifun(Dr_bands(1,2)*ones(size(CSRv)),CSRv) ); DR1_lo = DR1_lo(:);
DR1_hi = max( chifun(Dr_bands(1,1)*ones(size(CSRv)),CSRv), chifun(Dr_bands(1,2)*ones(size(CSRv)),CSRv) ); DR1_hi = DR1_hi(:);
DR2_lo = min( chifun(Dr_bands(2,1)*ones(size(CSRv)),CSRv), chifun(Dr_bands(2,2)*ones(size(CSRv)),CSRv) ); DR2_lo = DR2_lo(:);
DR2_hi = max( chifun(Dr_bands(2,1)*ones(size(CSRv)),CSRv), chifun(Dr_bands(2,2)*ones(size(CSRv)),CSRv) ); DR2_hi = DR2_hi(:);

% IC95 con clipping + medianas (LW=0.1)
hIC_B = gobjects(1,1);
for j = 1:numel(tests)
    [Clo,~,Chi] = ci_curve_bootstrap(CSRv, Dr_pts(j), bBoot, 'CSR', bOLS);
    Clo = Clo(:); Chi = Chi(:);

    if Dr_pts(j) >= Dr_bands(1,1) && Dr_pts(j) <= Dr_bands(1,2)
        Llo = DR1_lo; Lhi = DR1_hi;
    elseif Dr_pts(j) >= Dr_bands(2,1) && Dr_pts(j) <= Dr_bands(2,2)
        Llo = DR2_lo; Lhi = DR2_hi;
    else
        Llo = -inf(numel(CSRv),1); Lhi = inf(numel(CSRv),1);
    end

    Clo = max(Clo, Llo);
    Chi = min(Chi, Lhi);
    epsA = 1e-10; tight = Chi < Clo + epsA; Chi(tight) = Clo(tight) + epsA;

    h = fill_between(CSRv, Clo, Chi, cols(j,:), icAlpha, 'IC95 (todas las pruebas)');
    if j==1, hIC_B = h; else, set(h,'HandleVisibility','off'); end

    % mediana IC95 (LW=0.1)
    XcB   = [ones(numel(CSRv),1), log(Dr_pts(j))*ones(numel(CSRv),1), log(CSRv(:))];
    chiMedB = exp(prctile(XcB*bBoot.', 50, 2));
    plot(ax2, CSRv, chiMedB, 'k-', 'LineWidth', 0.1, 'HandleVisibility','off');
end

% Datos
pDataB = scatter(ax2, CSR_pts, chi_cal, 44, 'k','filled','DisplayName','Datos calibrados');
for i=1:numel(chi_cal)
    text(ax2, CSR_pts(i)+0.008*range(CSRv), chi_cal(i)+0.015*range(yl), labels(i), ...
         'FontName',fontName,'FontSize',7.5,'Clipping','on');
end

xlim(ax2,[min(CSRv) max(CSRv)]); ylim(ax2, yl)
xlabel(ax2,'$\mathrm{CSR}$ (-)','Interpreter','latex'); 
ylabel(ax2,'$\chi_{\max}$ (-)','Interpreter','latex')
legend(ax2,[hDR1 hDR2 hIC_B pDataB], ...
       {'$D_r$ 0.60–0.70 (ley)','$D_r$ 0.70–0.80 (ley)','IC95 (todas)','Datos'}, ...
       'Interpreter','latex','Location','northwest','Box','off');

%% ================== (C) LOO residual ==================
ax3 = nexttile(tlo,3); hold(ax3,'on'); box(ax3,'off')
xmin = min([chi_cal; pred_LOO])*0.95;  xmax = max([chi_cal; pred_LOO])*1.05;
e_rel = 100*(pred_LOO - chi_cal)./chi_cal;

pBand = patch(ax3, [xmin xmax xmax xmin], [-100*tol -100*tol 100*tol 100*tol], ...
              [0.93 0.55 0.55], 'FaceAlpha',0.35, 'EdgeColor','none', ...
              'DisplayName', sprintf('$\\pm %g\\,\\%%$ tolerancia (no CI)',100*tol));
yline(ax3,0,'k-','LineWidth',0.8,'HandleVisibility','off');
yline(ax3, 100*tol,'k--','LineWidth',0.8,'HandleVisibility','off');
yline(ax3,-100*tol,'k--','LineWidth',0.8,'HandleVisibility','off');

pDots = scatter(ax3, chi_cal, e_rel, 40, 'k','filled', 'DisplayName','$\hat{\chi}_{\max}^{\mathrm{LOO}}$');
for i=1:numel(chi_cal)
    dx = 0.012*diff([xmin xmax]);  dy = 0.012*max(5,range(e_rel));
    text(ax3, chi_cal(i)+dx, e_rel(i)+sign(e_rel(i))*abs(dy), sprintf('$\\mathrm{TCUI}\\,%s$', nums(i)), ...
         'Interpreter','latex','FontName',fontName,'FontSize',7.5,'Clipping','on');
end
legend(ax3,[pBand pDots], 'Interpreter','latex','Location','southeast','Box','off');
xlabel(ax3,'$\chi_{\max}$ calibrated (-)','Interpreter','latex')
ylabel(ax3,'Relative error $e_r$ (\%)','Interpreter','latex')
xlim(ax3,[xmin xmax]);
ylim(ax3,[min(min(e_rel)-5, -100*tol-5), max(max(e_rel)+5, 100*tol+5)]);
text(ax3, xmin, min(ylim(ax3))+0.04*diff(ylim(ax3)), ...
    sprintf('$e_r=100\\,(\\hat{\\chi}-\\chi)/\\chi$; bias=%.2f\\%%, RMSE_{rel}=%.2f\\%%, MAPE=%.2f\\%%', ...
            biasRel, RMSErel, MAPE), ...
    'Interpreter','latex','FontName',fontName,'FontSize',9);

%% (Opcional) imprimir IC95 puntuales en los puntos observados
fprintf('\nIC95%% puntuales en los puntos observados (χ_max):\n')
for j = 1:numel(tests)
    x0 = [1 log(Dr_pts(j)) log(CSR_pts(j))];
    [lo,md,hi] = ci_point_bootstrap(x0, bBoot);
    fprintf('%s: median=%.2f; IC95%%=[%.2f, %.2f]\n', upper(tests(j)), md, lo, hi);
end

%% -------------------- Sensibilidad a chi_cal(TCUI20) --------------------
% Si tus datos se llaman chi_exp, crea alias chi_cal
if exist('chi_cal','var')~=1 && exist('chi_exp','var')==1
    chi_cal = chi_exp;
end

% Índice del ensayo TCUI20
idx20 = find(strcmpi(tests,'tcui20'), 1);
if isempty(idx20)
    error('No se encontró TCUI20 en "tests".');
end

% Leverage local (NO usar 'h' para evitar choque con handles gráficos)
lev = diag( X / (X'*X) * X.' );   % vector n×1

% --- Estado base (con los valores actuales) ---
beta_base   = bOLS;           % [β0; β1; β2]
yhat_base   = X*beta_base;
res_base    = Y - yhat_base;
yhatLOO_b   = yhat_base - res_base./(1 - lev);
predLOO_b   = exp(yhatLOO_b);
e_rel_b     = 100*(predLOO_b - exp(Y))./exp(Y);
viol_base   = sum(abs(e_rel_b) > 100*tol);

% --- Cambio puntual en TCUI20: 31.5 -> 32.5 ---
chi_old = chi_cal(idx20);
chi_new = 32.5;                               % <— tu modificación
Delta   = log(chi_new) - log(chi_old);

XtX_inv = inv(X'*X);
beta_new = beta_base + XtX_inv * (X(idx20,:).') * Delta;

% Recalcular con Y modificado solo en TCUI20
Ynew       = Y;           Ynew(idx20) = log(chi_new);
yhat_new   = X*beta_new;
res_new    = Ynew - yhat_new;
yhatLOO_n  = yhat_new - res_new./(1 - lev);
predLOO_n  = exp(yhatLOO_n);
e_rel_n    = 100*(predLOO_n - exp(Ynew))./exp(Ynew);

MAPE_new    = mean(abs(e_rel_n));
RMSErel_new = sqrt(mean(e_rel_n.^2));
biasRel_new = mean(e_rel_n);
viol_new    = sum(abs(e_rel_n) > 100*tol);

fprintf('\n=== Sensibilidad puntual (TCUI20: %.1f -> %.1f) ===\n', chi_old, chi_new);
fprintf('Violaciones > ±%.0f%%%%: base=%d -> nuevo=%d\n', 100*tol, viol_base, viol_new);
fprintf('MAPE: %.2f%% -> %.2f%% | RMSE_rel: %.2f%% -> %.2f%% | Bias: %.2f%% -> %.2f%%\n', ...
        mean(abs(e_rel_b)), MAPE_new, sqrt(mean(e_rel_b.^2)), RMSErel_new, mean(e_rel_b), biasRel_new);

% (Opcional) ver el desplazamiento en parámetros (a,b,c)
a_base = exp(beta_base(1)); b_base = beta_base(2); c_base = beta_base(3);
a_new  = exp(beta_new(1));  b_new  = beta_new(2);  c_new  = beta_new(3);
fprintf('Parámetros: a=%.2f->%.2f | b=%.3f->%.3f | c=%.3f->%.3f\n', ...
        a_base, a_new, b_base, b_new, c_base, c_new);
end

%% ===================== Funciones locales =====================
function [Zlo,Zmd,Zhi] = ci_curve_bootstrap(Drv_or_CSRv, fixed, bBoot, mode, bOLS)
% Banda de CONFIANZA 95% (media en log) por bootstrap de parámetros (percentil).
    N = numel(Drv_or_CSRv);
    if strcmpi(mode,'Dr')
        x1 = log(Drv_or_CSRv(:)); x2 = log(fixed)*ones(N,1);
    elseif strcmpi(mode,'CSR')
        x1 = log(fixed)*ones(N,1); x2 = log(Drv_or_CSRv(:));
    else, error('mode: Dr o CSR'); end
    Xc = [ones(N,1) x1 x2];
    if ~isempty(bBoot)
        Yc  = Xc * bBoot.';                  % log de la MEDIA predicha
        CIy = prctile(Yc, [2.5 50 97.5], 2); % percentil
        Zlo = exp(CIy(:,1)); Zmd = exp(CIy(:,2)); Zhi = exp(CIy(:,3));
    else
        ymd = Xc*bOLS; Zmd = exp(ymd); Zlo = NaN(size(ymd)); Zhi = NaN(size(ymd));
    end
end

function [lo,md,hi] = ci_point_bootstrap(x0, bBoot)
% IC95 puntual (media en log) en x0, por bootstrap de parámetros.
    Yb = (x0*bBoot.').';             % (B x 1) log
    q  = prctile(Yb,[2.5 50 97.5]);
    lo = exp(q(1)); md = exp(q(2)); hi = exp(q(3));
end

function h = fill_between(x, ylo, yhi, col, alpha, legtxt)
% Sombreado robusto entre 2 curvas (sanea dimensiones y NaN).
    x   = x(:);  ylo = ylo(:);  yhi = yhi(:);
    if numel(x)~=numel(ylo) || numel(x)~=numel(yhi)
        error('fill_between: length mismatch: x=%d, ylo=%d, yhi=%d', ...
               numel(x), numel(ylo), numel(yhi));
    end
    bad = isnan(ylo) | isnan(yhi);
    if any(bad), x=x(~bad); ylo=ylo(~bad); yhi=yhi(~bad); end
    if isempty(x), h = gobjects(1); return; end
    h = patch([x; flipud(x)], [ylo; flipud(yhi)], col, ...
              'FaceAlpha', alpha, 'EdgeColor','none', 'DisplayName',legtxt, ...
              'HandleVisibility','on');  % permite proxy de leyenda
end
%%


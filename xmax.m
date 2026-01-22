%% ===== Calibración χ_max(Dr, CSR) con LOO + diagnóstico (Figura 1x3) =====
clear; close all; clc

% ---------- Datos ----------
tests   = ["tcui8","tcui9","tcui10","tcui20","tcui19","tcui18"]';
Dr_pts  = [0.62  0.68  0.61   0.79   0.78   0.78]';
CSR_pts = [0.20  0.25  0.30   0.50   0.40   0.30]';
chi_cal = [13.7  16.7  26.8   28.2   24.8   16.0]';   % χ_max (xmax)

% Parámetros de control
useRobust = false;      % true -> robustfit('huber') en log-espacio
tolPct    = 5;          % tolerancia para clasificar Over/Under (por ejemplo, 5 o 10)

% Etiquetas "TCUI n"
nums   = regexp(tests,'\d+','match','once');
labels = "TCUI " + string(nums);

% ---------- Ajuste en log-espacio: log(χ)=β0 + β1 log(Dr) + β2 log(CSR) ----------
Y = log(chi_cal);
X = [ones(size(Dr_pts)) log(Dr_pts) log(CSR_pts)];
if ~useRobust
    beta = X\Y;   % MCO
else
    [br,~] = robustfit([log(Dr_pts) log(CSR_pts)], Y, 'huber');  % br = [b0; b1; b2]
    beta = [br(1); br(2); br(3)];
end
a = exp(beta(1)); b = beta(2); c = beta(3);
chifun = @(Dr,CSR) a.*(Dr.^b).*(CSR.^c);

% ---------- Leave-One-Out coherente con el modo de ajuste ----------
n = numel(chi_cal);
pred_LOO = nan(n,1);
for i = 1:n
    m  = true(n,1); m(i) = false;
    if ~useRobust
        Xi = [ones(sum(m),1) log(Dr_pts(m)) log(CSR_pts(m))];
        Yi = log(chi_cal(m));
        bi = Xi\Yi;
    else
        [br_i,~] = robustfit([log(Dr_pts(m)) log(CSR_pts(m))], log(chi_cal(m)), 'huber');
        bi = [br_i(1); br_i(2); br_i(3)];
    end
    pred_LOO(i) = exp([1 log(Dr_pts(i)) log(CSR_pts(i))]*bi);
end

% ---------- Métricas en escala relativa ----------
rel      = (pred_LOO - chi_cal)./chi_cal;        % errores relativos adimensionales
MAPE     = 100*mean(abs(rel));
RMSE_rel = 100*sqrt(mean(rel.^2));
bias     = 100*mean(rel);
R2_id    = 1 - sum((pred_LOO-chi_cal).^2)/sum((chi_cal-mean(chi_cal)).^2);
P        = polyfit(chi_cal, pred_LOO, 1);        % paridad (pendiente/intercepto informativos)
m_lin    = P(1);  b_lin = P(2);

% ---------- Clasificación Over / Under / Within ±tol ----------
tol = tolPct/100;
class = strings(n,1);
class(rel >  tol) = "Over";
class(rel < -tol) = "Under";
class(abs(rel) <= tol) = sprintf("Within ±%d%%", tolPct);

% ---------- Reporte en consola (tabla ordenada por |error| desc) ----------
[~,ord] = sort(abs(rel),'descend');
T = table(tests(ord), Dr_pts(ord), CSR_pts(ord), chi_cal(ord), pred_LOO(ord), 100*rel(ord), class(ord), ...
    'VariableNames', {'Test','Dr','CSR','chi_cal','chi_LOO','e_rel_pct','Class'});
disp(T)
fprintf('chi_max fit -> a=%.4g, b=%.4g, c=%.4g | slope=%.2f, int=%.2f | R2_id=%.3f | MAPE=%.2f%% | RMSE_rel=%.2f%% | bias=%.2f%% | robust=%d\n', ...
        a,b,c,m_lin,b_lin,R2_id,MAPE,RMSE_rel,bias,useRobust);

% ---------- Rango de evaluación ----------
Drv  = linspace(0.55,0.85,400);
CSRv = linspace(0.20,0.50,400);

% ---------- Bandas por rangos (no son ICs) ----------
CSR_bands = [0.20 0.30; 0.30 0.50];     % panel vs Dr
Dr_bands  = [0.60 0.70; 0.70 0.80];     % panel vs CSR
colLow    = [0.20 0.45 0.90];           % azul
colHigh   = [0.85 0.30 0.35];           % rojo
alphaBand = 0.22;

% ---------- Estilo general ----------
set(0,'DefaultAxesFontName','Times New Roman','DefaultAxesFontSize',11);
set(0,'DefaultAxesTickDir','out','DefaultAxesLineWidth',0.5);
mk   = 44;  lwGrey = 0.9;

% ---------- Límites y rango χ ----------
[DRg, CSg] = ndgrid(Drv, CSRv);
CHI_all    = chifun(DRg, CSg);
yl = [min([CHI_all(:); chi_cal])*0.95, max([CHI_all(:); chi_cal])*1.05];

% ---------- Figura 1x3 ----------
fig = figure('Color','w','Units','centimeters','Position',[2 2 28 10]);
tlo = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');

% (1) χ_max vs Dr con bandas de CSR (bandas robustas con min/max)
ax1 = nexttile(tlo,1); hold(ax1,'on')
y1 = chifun(Drv, CSR_bands(1,1));  y2 = chifun(Drv, CSR_bands(1,2));
patch(ax1,[Drv, fliplr(Drv)], [min(y1,y2), fliplr(max(y1,y2))], colLow,  'FaceAlpha',alphaBand,'EdgeColor','none','DisplayName','0.20 \leq CSR < 0.30');
y1 = chifun(Drv, CSR_bands(2,1));  y2 = chifun(Drv, CSR_bands(2,2));
patch(ax1,[Drv, fliplr(Drv)], [min(y1,y2), fliplr(max(y1,y2))], colHigh, 'FaceAlpha',alphaBand,'EdgeColor','none','DisplayName','0.30 \leq CSR \leq 0.50');
for i=1:n   % iso-curvas grises por CSR de cada punto
    plot(ax1, Drv, chifun(Drv, CSR_pts(i)), 'Color',[0.6 0.6 0.6], 'LineWidth', lwGrey, 'HandleVisibility','off');
end
s1 = scatter(ax1, Dr_pts, chi_cal, mk, 'k','filled','DisplayName','Calibrated points');
dx = 0.006*(max(Drv)-min(Drv));  dy = 0.015*(yl(2)-yl(1));
for i=1:n
    text(Dr_pts(i)+dx, chi_cal(i)+dy, labels(i), 'FontName','Times New Roman','FontSize',7.5,'Clipping','on');
end
xlim(ax1,[min(Drv) max(Drv)]); ylim(ax1, yl)
xlabel(ax1,'$D_r$ (-)','Interpreter','latex'); ylabel(ax1,'$\chi_{\max}$ (-)','Interpreter','latex')
box(ax1,'off'); grid(ax1,'off')
legend(ax1,'show','Interpreter','latex','Location','northwest','NumColumns',1,'Box','off');

% (2) χ_max vs CSR con bandas de Dr (bandas robustas con min/max)
ax2 = nexttile(tlo,2); hold(ax2,'on')
y1 = chifun(Dr_bands(1,1)*ones(size(CSRv)), CSRv);
y2 = chifun(Dr_bands(1,2)*ones(size(CSRv)), CSRv);
patch(ax2,[CSRv, fliplr(CSRv)], [min(y1,y2), fliplr(max(y1,y2))], colLow,  'FaceAlpha',alphaBand,'EdgeColor','none','DisplayName','0.60 \leq D_r < 0.70');
y1 = chifun(Dr_bands(2,1)*ones(size(CSRv)), CSRv);
y2 = chifun(Dr_bands(2,2)*ones(size(CSRv)), CSRv);
patch(ax2,[CSRv, fliplr(CSRv)], [min(y1,y2), fliplr(max(y1,y2))], colHigh, 'FaceAlpha',alphaBand,'EdgeColor','none','DisplayName','0.70 \leq D_r \leq 0.80');
for i=1:n
    plot(ax2, CSRv, chifun(Dr_pts(i)*ones(size(CSRv)), CSRv), 'Color',[0.6 0.6 0.6], 'LineWidth', lwGrey, 'HandleVisibility','off');
end
s2 = scatter(ax2, CSR_pts, chi_cal, mk, 'k','filled','DisplayName','Calibrated points');
dx = 0.008*(max(CSRv)-min(CSRv));  dy = 0.015*(yl(2)-yl(1));
for i=1:n
    text(CSR_pts(i)+dx, chi_cal(i)+dy, labels(i), 'FontName','Times New Roman','FontSize',7.5,'Clipping','on');
end
xlim(ax2,[min(CSRv) max(CSRv)]); ylim(ax2, yl)
xlabel(ax2,'CSR (-)','Interpreter','latex'); ylabel(ax2,'$\chi_{\max}$ (-)','Interpreter','latex')
box(ax2,'off'); grid(ax2,'off')
legend(ax2,'show','Interpreter','latex','Location','northwest','NumColumns',1,'Box','off');

% (3) Paridad LOO con banda ±tolPct% y ajuste lineal
ax3 = nexttile(tlo,3); hold(ax3,'on')
xmin = min([chi_cal; pred_LOO])*0.95; xmax = max([chi_cal; pred_LOO])*1.05;
xb   = linspace(xmin, xmax, 300);
patch(ax3,[xb, fliplr(xb)], [(1 - tolPct/100)*xb, fliplr((1 + tolPct/100)*xb)], [0.75 0.90 0.80], 'FaceAlpha',0.30, 'EdgeColor','none','DisplayName',sprintf('\\pm %d\\%%',tolPct));
scatter(ax3, chi_cal, pred_LOO, mk*0.9, 'k', 'filled','DisplayName','LOO points');
% Etiquetas + clasificación color (Over rojo / Under azul / Within gris)
for i=1:n
    if rel(i) >  tol, col = [0.85 0.33 0.10];      % Over
    elseif rel(i) < -tol, col = [0.20 0.45 0.90];  % Under
    else, col = [0.40 0.40 0.40];                  % Within
    end
    plot(ax3, chi_cal(i), pred_LOO(i), 'o', 'MarkerSize', 7.5, 'MarkerFaceColor','none', 'MarkerEdgeColor', col, 'LineWidth',1.1, 'HandleVisibility','off');
    text(chi_cal(i)+0.012*(xmax-xmin), pred_LOO(i)+0.012*(xmax-xmin), labels(i), 'FontName','Times New Roman','FontSize',7.5,'Clipping','on');
end
% Línea de mejor ajuste (opcional informativa)
plot(ax3, [xmin xmax], polyval(P,[xmin xmax]), 'Color',[0.10 0.45 0.85], 'LineWidth',1.2, 'DisplayName',sprintf('$\\hat{y}=%.2f\\,x%+0.2f$', m_lin, b_lin));

axis(ax3,'equal'); xlim(ax3,[xmin xmax]); ylim(ax3,[xmin xmax])
xlabel(ax3,'$\chi_{\max}$ calibrated','Interpreter','latex')
ylabel(ax3,'$\hat{\chi}_{\max}^{\mathrm{LOO}}$','Interpreter','latex')
box(ax3,'off'); grid(ax3,'off')
legend(ax3,'show','Interpreter','latex', 'Location','southeast','Box','off');




%% ===== Calibración z_max(D_r, CSR) — Figura 1×3 con residuo relativo =====
% Ajuste en log-espacio; métricas relativas; panel 3 = error relativo (%)

clear; close all; clc

% -------------------- Parámetros --------------------
useRobust  = false;          % true -> robustfit (Huber); false -> OLS
tol        = 0.10;           % banda de tolerancia del panel 3 (±10%)
fontName   = 'Times New Roman';
fontSize   = 11;

% -------------------- Datos --------------------
tests  = ["tcui8","tcui9","tcui10","tcui20","tcui19","tcui18"]';
Dr_pts = [0.62  0.68  0.61  0.79  0.78  0.78]';
CSR_pts= [0.20  0.25  0.30  0.50  0.40  0.30]';
% z_cal  = [7.0   9.2   18.3  28.0  19.5   9.3 ]';   % zmax manual
z_cal  = [7.1   9.2   18.5  28.2  19.6   9.8 ]';   % zmax manual
% Etiquetas "TCUI n"
nums   = regexp(tests,'\d+','match','once');
labels = "TCUI " + string(nums);

% -------------------- Ajuste (log-espacio) --------------------
% log(z) = β0 + β1*log(Dr) + β2*log(CSR)  <=>  z = a*Dr^b*CSR^c
Y = log(z_cal);
X = [ones(size(Dr_pts))  log(Dr_pts)  log(CSR_pts)];

if ~useRobust
    beta = X\Y;                                   % OLS
else
    [br,~] = robustfit([log(Dr_pts) log(CSR_pts)], Y, 'huber');
    beta   = [br(1); br(2); br(3)];               % intercepto + coeficientes
end

a = exp(beta(1)); b = beta(2); c = beta(3);
zfun = @(Dr,CSR) a.*(Dr.^b).*(CSR.^c);

% -------------------- Leave-One-Out (LOO) --------------------
pred_LOO = nan(size(z_cal));
for i = 1:numel(z_cal)
    msk = true(size(z_cal)); msk(i) = false;
    Xi = [ones(sum(msk),1) log(Dr_pts(msk)) log(CSR_pts(msk))];
    Yi = log(z_cal(msk));
    if ~useRobust
        bi = Xi\Yi;
    else
        br_i = robustfit([log(Dr_pts(msk)) log(CSR_pts(msk))], Yi, 'huber');
        bi   = [br_i(1); br_i(2); br_i(3)];
    end
    pred_LOO(i) = exp([1 log(Dr_pts(i)) log(CSR_pts(i))]*bi);
end

% -------------------- Métricas relativas (consistentes con log-OLS) ------
rel     = (pred_LOO - z_cal)./z_cal;
MAPE    = 100*mean(abs(rel));
RMSErel = 100*sqrt(mean(rel.^2));
biasRel = 100*mean(rel);
R2_id   = 1 - sum((pred_LOO - z_cal).^2)/sum((z_cal - mean(z_cal)).^2);

% (para referencia: recta libre del parity, NO necesaria en residual)
P     = polyfit(z_cal, pred_LOO, 1);
slope = P(1); interc = P(2);

fprintf('zmax: a=%.4g b=%.4g c=%.4g | slope=%.2f, int=%.2f | R2_id=%.3f | MAPE=%.2f%%, RMSErel=%.2f%%, bias=%.2f%%\n',...
        a,b,c,slope,interc,R2_id,MAPE,RMSErel,biasRel);

% -------------------- Rango de evaluación --------------------
Drv  = linspace(0.55,0.85,400);
CSRv = linspace(0.20,0.50,400);

% Bandas (no ICs)
CSR_bands = [0.20 0.30; 0.30 0.50];     % panel vs Dr
Dr_bands  = [0.60 0.70; 0.70 0.80];     % panel vs CSR

% Paleta
colLow    = [0.00 0.62 0.45];   % verde azulado
colHigh   = [0.88 0.42 0.38];   % coral suave
alphaBand = 0.20;

% -------------------- Estilo global --------------------
set(0,'DefaultAxesFontName',fontName,'DefaultAxesFontSize',fontSize);
set(0,'DefaultAxesTickDir','out','DefaultAxesLineWidth',0.5);
mk = 44;  lwGrey = 0.9;

% -------------------- Límites comunes (y) --------------------
[DRg, CSg] = ndgrid(Drv, CSRv);
Z_all      = zfun(DRg, CSg);
yl = [min([Z_all(:); z_cal])*0.95, max([Z_all(:); z_cal])*1.05];

% -------------------- Figura 1×3 --------------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 28 10]);
tlo = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');

% (1) zmax vs Dr con bandas de CSR (rango, no ICs)
ax1 = nexttile(tlo,1); hold(ax1,'on')

y1 = zfun(Drv, CSR_bands(1,1)); y2 = zfun(Drv, CSR_bands(1,2));
patch(ax1,[Drv, fliplr(Drv)], [min(y1,y2), fliplr(max(y1,y2))], colLow,  ...
      'FaceAlpha',alphaBand,'EdgeColor','none','DisplayName','0.20 \leq CSR < 0.30');
y1 = zfun(Drv, CSR_bands(2,1)); y2 = zfun(Drv, CSR_bands(2,2));
patch(ax1,[Drv, fliplr(Drv)], [min(y1,y2), fliplr(max(y1,y2))], colHigh, ...
      'FaceAlpha',alphaBand,'EdgeColor','none','DisplayName','0.30 \leq CSR \leq 0.50');

for i=1:numel(CSR_pts)
    plot(ax1, Drv, zfun(Drv, CSR_pts(i)), 'Color',[0.6 0.6 0.6], ...
         'LineWidth', lwGrey, 'HandleVisibility','off');
end

scatter(ax1, Dr_pts, z_cal, mk, 'k','filled','DisplayName','Calibrated points');
dx = 0.006*range(Drv);  dy = 0.015*range(yl);
for i=1:numel(z_cal)
    text(ax1, Dr_pts(i)+dx, z_cal(i)+dy, labels(i), ...
         'FontName',fontName,'FontSize',7.5,'Clipping','on');
end
xlim(ax1,[min(Drv) max(Drv)]); ylim(ax1, yl)
xlabel(ax1,'D_r (-)','Interpreter','latex'); ylabel(ax1,'$z_{\max}$ (-)','Interpreter','latex')
box(ax1,'off'); grid(ax1,'off')
legend(ax1,'show','Interpreter','latex','Location','northwest','Box','off');

% (2) zmax vs CSR con bandas de Dr (rango, no ICs)
ax2 = nexttile(tlo,2); hold(ax2,'on')

y1 = zfun(Dr_bands(1,1)*ones(size(CSRv)), CSRv);
y2 = zfun(Dr_bands(1,2)*ones(size(CSRv)), CSRv);
patch(ax2,[CSRv, fliplr(CSRv)], [min(y1,y2), fliplr(max(y1,y2))], colLow, ...
      'FaceAlpha',alphaBand,'EdgeColor','none','DisplayName','0.60 \leq D_r < 0.70');
y1 = zfun(Dr_bands(2,1)*ones(size(CSRv)), CSRv);
y2 = zfun(Dr_bands(2,2)*ones(size(CSRv)), CSRv);
patch(ax2,[CSRv, fliplr(CSRv)], [min(y1,y2), fliplr(max(y1,y2))], colHigh,...
      'FaceAlpha',alphaBand,'EdgeColor','none','DisplayName','0.70 \leq D_r \leq 0.80');

for i=1:numel(Dr_pts)
    plot(ax2, CSRv, zfun(Dr_pts(i)*ones(size(CSRv)), CSRv), 'Color',[0.6 0.6 0.6], ...
         'LineWidth', lwGrey, 'HandleVisibility','off');
end

scatter(ax2, CSR_pts, z_cal, mk, 'k','filled','DisplayName','Calibrated points');
dx = 0.008*range(CSRv);  dy = 0.015*range(yl);
for i=1:numel(z_cal)
    text(ax2, CSR_pts(i)+dx, z_cal(i)+dy, labels(i), ...
         'FontName',fontName,'FontSize',7.5,'Clipping','on');
end
xlim(ax2,[min(CSRv) max(CSRv)]); ylim(ax2, yl)
xlabel(ax2,'CSR (-)','Interpreter','latex'); ylabel(ax2,'$z_{\max}$ (-)','Interpreter','latex')
box(ax2,'off'); grid(ax2,'off')
legend(ax2,'show','Interpreter','latex','Location','northwest','Box','off');

% (3) Residuo relativo (%)
ax3 = nexttile(tlo,3); hold(ax3,'on')
xmin = min([z_cal; pred_LOO])*0.95;
xmax = max([z_cal; pred_LOO])*1.05;

e_rel = 100*(pred_LOO - z_cal)./z_cal;

pBand = patch(ax3, [xmin xmax xmax xmin], [-100*tol -100*tol 100*tol 100*tol], ...
              [0.75 0.90 0.80], 'FaceAlpha',0.35, 'EdgeColor','none');
yline(ax3,0,'k-','LineWidth',0.8,'HandleVisibility','off');
yline(ax3, 100*tol,'k--','LineWidth',0.8,'HandleVisibility','off');
yline(ax3,-100*tol,'k--','LineWidth',0.8,'HandleVisibility','off');

pDots = scatter(ax3, z_cal, e_rel, mk*0.9, 'k','filled');

% Etiquetas y leyenda rigurosa
xr = diff([xmin xmax]); yr = diff([min(e_rel) max(e_rel)]);
dx3 = 0.012*xr; dy3 = 0.012*max(5,abs(yr));
for i=1:numel(z_cal)
    text(ax3, z_cal(i)+dx3, e_rel(i)+sign(e_rel(i))*abs(dy3), ...
        sprintf('$\\mathrm{TCUI}\\,%s$', nums(i)), ...
        'Interpreter','latex','FontName',fontName,'FontSize',7.5,'Clipping','on');
end
legend(ax3,[pBand pDots], {sprintf('$\\pm %g\\,\\%%$ tolerance (not CI)',100*tol), ...
       '$\\hat{z}_{\\max}^{\\mathrm{LOO}}$'}, 'Interpreter','latex', ...
       'Location','southeast','Box','off');

% Texto breve (definición de e_r y métricas)
yr = ylim(ax3); xr = xlim(ax3);
text(ax3, xr(1), yr(1)+0.04*diff(yr), ...
    sprintf('$e_r=100\\,(\\hat z-z)/z$; bias=%.2f\\%%, RMSE_{rel}=%.2f\\%%, MAPE=%.2f\\%%', ...
            biasRel, RMSErel, MAPE), ...
    'Interpreter','latex','FontName',fontName,'FontSize',9);

xlim(ax3,[xmin xmax]);
ylim(ax3,[min(min(e_rel)-5, -100*tol-5), max(max(e_rel)+5, 100*tol+5)]);
xlabel(ax3,'$z_{\\max}$ calibrated (-)','Interpreter','latex')
ylabel(ax3,'Relative error $e_r$ (\%)','Interpreter','latex')
box(ax3,'off'); grid(ax3,'off')

% -------------------- Export (opcional) --------------------
% exportgraphics(fig,'zmax_calibration_residual.pdf','Resolution',600);

%%

%% ===== zmax_uncertainty_all_IC95.m =====
% Ley: log(z) = β0 + β1 log(D_r) + β2 log(CSR)
% (a) z vs D_r: medianas por CSR + IC95% (todas las pruebas)
% (b) z vs CSR: medianas por D_r + IC95% (todas las pruebas) con clipping a rangos de D_r
% (c) error relativo LOO con banda operativa ±5 %

clear; close all; clc; rng(1)

%% -------------------- Estilo / opciones --------------------
tol        = 0.05;                       % banda operativa (panel c)
fontName   = 'Times New Roman'; fontSize = 11;
set(0,'DefaultAxesFontName',fontName,'DefaultAxesFontSize',fontSize);
set(0,'DefaultAxesTickDir','out','DefaultAxesLineWidth',0.5);
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultAxesTickLabelInterpreter','latex');

icAlpha    = 0.20;                        % transparencia IC95 (bandas)
medLW      = 0.10;                        % grosor de las medianas
colMustard = [0.89 0.67 0.16];            % franjas-ley CSR / D_r
colGreen   = [0.00 0.62 0.45];
rangeAlpha = 0.14;

% Rangos operativos de referencia (NO IC)
CSR_bands = [0.20 0.30; 0.30 0.50];       % para panel (A)
Dr_bands  = [0.60 0.70; 0.70 0.80];       % para panel (B)

%% -------------------- Datos (corregidos) --------------------
tests   = ["tcui8","tcui9","tcui10","tcui20","tcui19","tcui18"]';
Dr_pts  = [0.62  0.68  0.61  0.79  0.78  0.78]';
CSR_pts = [0.20  0.25  0.30  0.50  0.40  0.30]';
z_cal   = [6.8  9.1   18.7  34.8  19.4   10]';
% z_cal   = [6.8  9.1   18.7  25.3  19.4   10]';
nums    = regexp(tests,'\d+','match','once');
labels  = "TCUI " + string(nums);

% Diseño en log
Y = log(z_cal);
X = [ones(size(Dr_pts))  log(Dr_pts)  log(CSR_pts)];
[n,p] = size(X);                          % n=6, p=3

%% -------------------- OLS, LOO exacto, métricas --------------------
bOLS  = (X'*X)\(X'*Y);                    % [β0; β1; β2]
yhat  = X*bOLS;
res   = Y - yhat;
H     = X/(X'*X)*X';  h = diag(H);

a = exp(bOLS(1)); b = bOLS(2); c = bOLS(3);
zfun = @(Dr,CSR) a.*(Dr.^b).*(CSR.^c);

% LOO exacto (panel c)
yhat_LOO = Y - res./(1 - h);
pred_LOO = exp(yhat_LOO);
rel      = (pred_LOO - z_cal)./z_cal;
MAPE     = 100*mean(abs(rel));
RMSErel  = 100*sqrt(mean(rel.^2));
biasRel  = 100*mean(rel);

%% -------------------- Bootstrap de parámetros (wild Rademacher) --------------------
B = 5000;                                  % 2000–5000 ok para n=6
u_tilde = res ./ (1 - h);                  % residuo corregido
bBoot   = zeros(B,3);
XtX_inv = inv(X'*X);                       % eficiencia
for k = 1:B
    v    = sign(rand(n,1)-0.5);            % ±1 (Rademacher)
    Yb   = X*bOLS + u_tilde .* v;          % bootstrap en log
    bBoot(k,:) = XtX_inv*(X'*Yb);          % re-ajuste OLS
end

%% -------------------- Rangos y límites --------------------
Drv  = linspace(0.55,0.85,400);
CSRv = linspace(0.20,0.50,400);
[DRg, CSg] = ndgrid(Drv, CSRv); %#ok<ASGLU>
Z_all = zfun(DRg, CSg);
yl = [min([Z_all(:); z_cal])*0.95, max([Z_all(:); z_cal])*1.05];
dy = diff(yl);

%% -------------------- FIGURA 1×3 --------------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 28 10]);
tlo = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');

%% ========== (A) zmax vs D_r ==========
ax1 = nexttile(tlo,1); hold(ax1,'on'); box(ax1,'off')

% Franjas de CSR (ley, no IC) — al fondo
y1 = zfun(Drv, CSR_bands(1,1)); y2 = zfun(Drv, CSR_bands(1,2));
hCSR1 = fill_between(Drv, min(y1,y2), max(y1,y2), colMustard, rangeAlpha, 'CSR 0.20–0.30 (law)');
y1 = zfun(Drv, CSR_bands(2,1)); y2 = zfun(Drv, CSR_bands(2,2));
hCSR2 = fill_between(Drv, min(y1,y2), max(y1,y2), colGreen,   rangeAlpha, 'CSR 0.30–0.50 (law)');

uistack([hCSR1 hCSR2],'bottom');

% IC95 de todos los ensayos (escala de grises)
grayLevels = linspace(0.85, 0.30, numel(tests));
cols = repmat(grayLevels', 1, 3);
hIC  = gobjects(1,1); % proxy

for j = 1:numel(tests)
    [Zlo,~,Zhi] = ci_curve_bootstrap(Drv, CSR_pts(j), bBoot, 'Dr', bOLS);
    h = fill_between(Drv, Zlo, Zhi, cols(j,:), icAlpha, 'IC95\% (all tests)');
    if j==1, hIC = h; else, set(h,'HandleVisibility','off'); end
end

% Medianas (todas) con grosor 0.1
for j = 1:numel(tests)
    XcA   = [ones(numel(Drv),1), log(Drv(:)), log(CSR_pts(j))*ones(numel(Drv),1)];
    zMedA = exp(prctile(XcA*bBoot.', 50, 2));
    plot(ax1, Drv, zMedA, 'k-', 'LineWidth', medLW, 'HandleVisibility','off');
end

% Datos
pData = scatter(ax1, Dr_pts, z_cal, 44, 'k','filled','DisplayName','\mathrm{Data}');
for i=1:numel(z_cal)
    text(ax1, Dr_pts(i)+0.006*(Drv(end)-Drv(1)), z_cal(i)+0.015*dy, labels(i), ...
         'FontName',fontName,'FontSize',7.5,'Clipping','on');
end
xlim(ax1,[min(Drv) max(Drv)]); ylim(ax1, yl)
xlabel(ax1,'$D_r$ (-)'); 
ylabel(ax1,'$z_{\max}$ (-)')
legend(ax1,[hCSR1 hCSR2 hIC pData], ...
       {'$\mathrm{CSR}\in[0.20,\,0.30]$ (law)', ...
        '$\mathrm{CSR}\in[0.30,\,0.50]$ (law)', ...
        'IC95\,\% (all tests)', ...
        '\mathrm{Data}'}, ...
       'Location','northwest','Box','off');

%% ========== (B) zmax vs CSR (con clipping a rangos de D_r) ==========
ax2 = nexttile(tlo,2); hold(ax2,'on'); box(ax2,'off')

% Franjas de D_r (ley, no IC) — al fondo
y1 = zfun(Dr_bands(1,1)*ones(size(CSRv)), CSRv);
y2 = zfun(Dr_bands(1,2)*ones(size(CSRv)), CSRv);
hDR1 = fill_between(CSRv, min(y1,y2), max(y1,y2), colMustard, rangeAlpha, '$D_r\in[0.60,\,0.70]$ (law)');
y1 = zfun(Dr_bands(2,1)*ones(size(CSRv)), CSRv);
y2 = zfun(Dr_bands(2,2)*ones(size(CSRv)), CSRv);
hDR2 = fill_between(CSRv, min(y1,y2), max(y1,y2), colGreen,   rangeAlpha, '$D_r\in[0.70,\,0.80]$ (law)');
uistack([hDR1 hDR2],'bottom');

% Pre-cómputo de límites por franja (columnas)
DR1_lo = min( zfun(Dr_bands(1,1)*ones(size(CSRv)),CSRv), zfun(Dr_bands(1,2)*ones(size(CSRv)),CSRv) ); DR1_lo = DR1_lo(:);
DR1_hi = max( zfun(Dr_bands(1,1)*ones(size(CSRv)),CSRv), zfun(Dr_bands(1,2)*ones(size(CSRv)),CSRv) ); DR1_hi = DR1_hi(:);
DR2_lo = min( zfun(Dr_bands(2,1)*ones(size(CSRv)),CSRv), zfun(Dr_bands(2,2)*ones(size(CSRv)),CSRv) ); DR2_lo = DR2_lo(:);
DR2_hi = max( zfun(Dr_bands(2,1)*ones(size(CSRv)),CSRv), zfun(Dr_bands(2,2)*ones(size(CSRv)),CSRv) ); DR2_hi = DR2_hi(:);

% IC95 con clipping según D_r del ensayo
hIC2 = gobjects(1,1); % proxy
for j = 1:numel(tests)
    [Zlo,~,Zhi] = ci_curve_bootstrap(CSRv, Dr_pts(j), bBoot, 'CSR', bOLS);
    Zlo = Zlo(:);  Zhi = Zhi(:);

    if Dr_pts(j) >= Dr_bands(1,1) && Dr_pts(j) <= Dr_bands(1,2)
        Llo = DR1_lo;  Lhi = DR1_hi;
    elseif Dr_pts(j) >= Dr_bands(2,1) && Dr_pts(j) <= Dr_bands(2,2)
        Llo = DR2_lo;  Lhi = DR2_hi;
    else
        Llo = -inf(numel(CSRv),1);  Lhi =  inf(numel(CSRv),1);
    end

    ZloC = max(Zlo, Llo);
    ZhiC = min(Zhi, Lhi);

    epsA = 1e-10; tight = ZhiC < ZloC + epsA;
    ZhiC(tight) = ZloC(tight) + epsA;

    h = fill_between(CSRv, ZloC, ZhiC, cols(j,:), icAlpha, 'IC95\,\% (all tests)');
    if j==1, hIC2 = h; else, set(h,'HandleVisibility','off'); end
end

% Medianas (todas) con grosor 0.1
for j = 1:numel(tests)
    XcB   = [ones(numel(CSRv),1), log(Dr_pts(j))*ones(numel(CSRv),1), log(CSRv(:))];
    zMedB = exp(prctile(XcB*bBoot.', 50, 2));
    plot(ax2, CSRv, zMedB, 'k-', 'LineWidth', medLW, 'HandleVisibility','off');
end

% Datos
pData2 = scatter(ax2, CSR_pts, z_cal, 44, 'k','filled','DisplayName','\mathrm{Data}');
for i=1:numel(z_cal)
    text(ax2, CSR_pts(i)+0.008*(CSRv(end)-CSRv(1)), z_cal(i)+0.015*dy, labels(i), ...
         'FontName',fontName,'FontSize',7.5,'Clipping','on');
end
xlim(ax2,[min(CSRv) max(CSRv)]); ylim(ax2, yl)
xlabel(ax2,'$\mathrm{CSR}$ (-)'); 
ylabel(ax2,'$z_{\max}$ (-)')
legend(ax2,[hDR1 hDR2 hIC2 pData2], ...
       {'$D_r\in[0.60,\,0.70]$ (law)', ...
        '$D_r\in[0.70,\,0.80]$ (law)', ...
        'IC95\,\% (all tests)', ...
        '\mathrm{Data}'}, ...
       'Location','northwest','Box','off');

%% ========== (C) Residuo relativo ==========
ax3 = nexttile(tlo,3); hold(ax3,'on'); box(ax3,'off')
xmin = min([z_cal; pred_LOO])*0.95;  xmax = max([z_cal; pred_LOO])*1.05;
e_rel = 100*(pred_LOO - z_cal)./z_cal;
pBand = patch(ax3, [xmin xmax xmax xmin], [-100*tol -100*tol 100*tol 100*tol], ...
              [0.93 0.55 0.55], 'FaceAlpha',0.35, 'EdgeColor','none');
yline(ax3,0,'k-','LineWidth',0.8,'HandleVisibility','off');
yline(ax3, 100*tol,'k--','LineWidth',0.8,'HandleVisibility','off');
yline(ax3,-100*tol,'k--','LineWidth',0.8,'HandleVisibility','off');
pDots = scatter(ax3, z_cal, e_rel, 40, 'k','filled');
for i=1:numel(z_cal)
    dx = 0.012*(xmax-xmin);  dy_txt = 0.012*max(5,range(e_rel));
    text(ax3, z_cal(i)+dx, e_rel(i)+sign(e_rel(i))*abs(dy_txt), sprintf('$\\mathrm{TCUI}\\,%s$', nums(i)), ...
         'FontName',fontName,'FontSize',7.5,'Clipping','on');
end
legend(ax3,[pBand pDots], ...
       {sprintf('$\\pm %g\\,\\%%$ tolerance',100*tol), '$\\hat{z}_{\\max}^{\\mathrm{LOO}}$'}, ...
       'Location','southeast','Box','off');
xlabel(ax3,'$z_{\\max}$ calibrated (-)')
ylabel(ax3,'Relative error $e_r$ (\%)')
xlim(ax3,[xmin xmax]);
ylim(ax3,[min(min(e_rel)-5, -100*tol-5), max(max(e_rel)+5, 100*tol+5)]);

%% (Opcional) IC95 puntuales en los puntos observados
fprintf('\nIC95%% puntuales en los puntos observados:\n')
for j = 1:numel(tests)
    x0 = [1 log(Dr_pts(j)) log(CSR_pts(j))];
    [lo,md,hi] = ci_point_bootstrap(x0, bBoot);
    fprintf('%s: median=%.2f; IC95%%=[%.2f, %.2f]\n', upper(tests(j)), md, lo, hi);
end

%% ===================== Funciones locales =====================
function [Zlo,Zmd,Zhi] = ci_curve_bootstrap(Drv_or_CSRv, fixed, bBoot, mode, bOLS)
% Banda de CONFIANZA 95% (media en log) por bootstrap de parámetros (percentil).
    N = numel(Drv_or_CSRv);
    if strcmpi(mode,'Dr')
        x1 = log(Drv_or_CSRv(:)); x2 = log(fixed)*ones(N,1);
    elseif strcmpi(mode,'CSR')
        x1 = log(fixed)*ones(N,1); x2 = log(Drv_or_CSRv(:));
    else
        error('mode: Dr o CSR');
    end
    Xc = [ones(N,1) x1 x2];
    if ~isempty(bBoot)
        Yc  = Xc * bBoot.';                  % log de la MEDIA predicha (bootstrap)
        CIy = prctile(Yc, [2.5 50 97.5], 2); % IC95% (percentil)
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
% Sombreado entre 2 curvas con saneo de dimensiones/NaN
    x   = x(:); ylo = ylo(:); yhi = yhi(:);
    if numel(x)~=numel(ylo) || numel(x)~=numel(yhi)
        error('fill_between: length mismatch: x=%d, ylo=%d, yhi=%d', ...
               numel(x), numel(ylo), numel(yhi));
    end
    bad = isnan(ylo) | isnan(yhi);
    if any(bad), x=x(~bad); ylo=ylo(~bad); yhi=yhi(~bad); end
    if isempty(x), h = gobjects(1); return; end
    h = patch([x; flipud(x)], [ylo; flipud(yhi)], col, ...
              'FaceAlpha', alpha, 'EdgeColor','none', 'DisplayName',legtxt);
end

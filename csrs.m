%% CSR–Nini (Fig.15-like) — 3x2 legend (symbols only), LW=0.1 everywhere, vertical grid + CSR horizontals
clear; clc;

% -------------------- Datos (usuario) --------------------
% bb (Bio Bío)
% CSR_exp_bb = [0.2, 0.2, 0.33, 0.17];%qamp/2p0
CSR_exp_bb = [0.4, 0.4, 0.66, 0.34];%qamp/p0
N_exp_bb   = [23,   21,   4,    29];

CSR_sim_bb = [0.4, 0.4, 0.66, 0.34];
N_sim_bb   = [5,    20,   3,    28];

% ks60 (Karlsruhe ~60%)
% CSR_exp_ks60 = [0.1, 0.125, 0.15];%qamp/2p0
CSR_exp_ks60 = [0.2, 0.25, 0.3];%qamp/p0
N_exp_ks60   = [249, 100,  15];

CSR_sim_ks60 = [0.2, 0.25, 0.3];
N_sim_ks60   = [251, 80,   16];

% ks70 (Karlsruhe ~70–80%)
% CSR_exp_ks70 = [0.15, 0.2, 0.25];%qamp/
CSR_exp_ks70 = [0.3 0.4 0.5];%qamp/2p0
N_exp_ks70   = [54,  15,  6];

CSR_sim_ks70 = [0.3 0.4 0.5];
N_sim_ks70   = [64,  14,  8];

% -------------------- Parámetros globales de estilo --------------------
LW = 0.1;  % <<<<<< TODAS las líneas con grosor 0.1

% -------------------- Figura / estilo --------------------
figure('Color','w'); ax = axes; hold(ax,'on');

set(ax,'XScale','log', ...
       'Box','on', ...
       'TickDir','out', ...
       'LineWidth',LW, ...
       'FontName','Helvetica', ...
       'FontSize',11, ...
       'Layer','top', ...
       'GridLineStyle',':', ...
       'MinorGridLineStyle',':');

xlabel(ax,'Number of cycles, $N_{f}$','Interpreter','latex');
ylabel(ax,'Amplitude-pressure ratio, $q_{amp}/(p_0)$','Interpreter','latex');

xlim(ax,[1 1000]);
ylim(ax,[0.05 0.7]);

% -------------------- Retícula vertical (X): punteada + minor --------------------
ax.XGrid      = 'off';
ax.XMinorGrid = 'off';
ax.XMinorTick = 'off';

% Horizontal base apagada (las horizontales las agregamos por CSR)
ax.YGrid      = 'off';
ax.YMinorGrid = 'off';

ax.GridAlpha      = 0.20;
ax.MinorGridAlpha = 0.12;

% Líneas verticales continuas en x=10 y x=100 (LW=0.1)
xl10  = xline(ax,10 ,'--','LineWidth',LW,'Color',[0.15 0.15 0.15]);  xl10.HandleVisibility  = 'off';
xl100 = xline(ax,100,'--','LineWidth',LW,'Color',[0.15 0.15 0.15]);  xl100.HandleVisibility = 'off';

% -------------------- Líneas horizontales en los CSR dateados --------------------
CSR_levels = unique([CSR_exp_bb(:);CSR_sim_bb(:); ...
                     CSR_exp_ks60(:);CSR_sim_ks60(:); ...
                     CSR_exp_ks70(:);CSR_sim_ks70(:)], 'sorted');

for k = 1:numel(CSR_levels)
    yl = yline(ax, CSR_levels(k), ':', 'LineWidth', LW, 'Color', [0.75 0.75 0.75]);
    yl.HandleVisibility = 'off';
end

% -------------------- Colores por arena --------------------
cBB = [0 0 0];                 % Bio Bío = negro
cKS = [0 0.4470 0.7410];        % Karlsruhe = azul

% -------------------- Marcadores por grupo --------------------
mkBB   = '^';   % BBS
mkKS60 = 'o';   % KFS Dr = 60–68%
mkKS70 = 's';   % KFS Dr = 78–79%

% Helper
plotCurve = @(N,CSR,mk,ls,col,faceCol) localPlotCurve(ax,N,CSR,mk,ls,col,faceCol,LW);

% -------------------- Curvas reales --------------------
% Experimental: vacío + sólido
plotCurve(N_exp_bb,   CSR_exp_bb,   mkBB,   '-',  cBB, 'none');
plotCurve(N_exp_ks60, CSR_exp_ks60, mkKS60, '-',  cKS, 'none');
plotCurve(N_exp_ks70, CSR_exp_ks70, mkKS70, '-',  cKS, 'none');

% Model: relleno + discontinuo
plotCurve(N_sim_bb,   CSR_sim_bb,   mkBB,   '--', cBB, cBB);
plotCurve(N_sim_ks60, CSR_sim_ks60, mkKS60, '--', cKS, cKS);
plotCurve(N_sim_ks70, CSR_sim_ks70, mkKS70, '--', cKS, cKS);

% -------------------- Leyenda 3x2: sólo símbolos (sin texto) --------------------
% Columna izquierda  = Experimental (vacío, sólido)
% Columna derecha    = Model        (relleno, discontinuo)

hE_ks60 = plot(ax,nan,nan,'-','Color',cKS,'LineWidth',LW, ...
    'Marker',mkKS60,'MarkerSize',7,'MarkerFaceColor','none','MarkerEdgeColor',cKS);

hE_ks70 = plot(ax,nan,nan,'-','Color',cKS,'LineWidth',LW, ...
    'Marker',mkKS70,'MarkerSize',7,'MarkerFaceColor','none','MarkerEdgeColor',cKS);

hE_bb   = plot(ax,nan,nan,'-','Color',cBB,'LineWidth',LW, ...
    'Marker',mkBB,'MarkerSize',7,'MarkerFaceColor','none','MarkerEdgeColor',cBB);

hM_ks60 = plot(ax,nan,nan,'--','Color',cKS,'LineWidth',LW, ...
    'Marker',mkKS60,'MarkerSize',7,'MarkerFaceColor',cKS,'MarkerEdgeColor',cKS);

hM_ks70 = plot(ax,nan,nan,'--','Color',cKS,'LineWidth',LW, ...
    'Marker',mkKS70,'MarkerSize',7,'MarkerFaceColor',cKS,'MarkerEdgeColor',cKS);

hM_bb   = plot(ax,nan,nan,'--','Color',cBB,'LineWidth',LW, ...
    'Marker',mkBB,'MarkerSize',7,'MarkerFaceColor',cBB,'MarkerEdgeColor',cBB);

% Con NumColumns=2, MATLAB llena column-major:
% [Exp fila1, Exp fila2, Exp fila3, Model fila1, Model fila2, Model fila3]
lgd = legend([hE_ks60,hE_ks70,hE_bb,  hM_ks60,hM_ks70,hM_bb], ...
    {'$\,$','$\,$','$\,$',  '$\,$','$\,$','$\,$'}, ...
    'NumColumns',2,'Location','northeast','Box','on','Interpreter','latex');

lgd.LineWidth = LW;
lgd.EdgeColor = [0 0 0];
lgd.Color     = [1 1 1];
lgd.ItemTokenSize = [14 9];

% -------------------- Función local --------------------
function h = localPlotCurve(ax,N,CSR,marker,ls,edgeCol,faceCol,LW)
    [N2,ix] = sort(N(:));
    CSR2 = CSR(ix);

    h = plot(ax,N2,CSR2, ...
        'LineStyle',ls,'LineWidth',LW,'Color',edgeCol, ...
        'Marker',marker,'MarkerSize',6, ...
        'MarkerEdgeColor',edgeCol, ...
        'MarkerFaceColor',faceCol);

    h.HandleVisibility = 'off'; % leyenda sólo con dummies
end



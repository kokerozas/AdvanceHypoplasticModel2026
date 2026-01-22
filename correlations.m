function correlations
% Ajustes vs D_r con IC95% (sin PI) + diagnósticos de residuales

% -------------------- Datos --------------------
TCUI     = {'tcui8','tcui9','tcui10','tcui20','tcui19','tcui18'};
Dr       = [0.62 0.68 0.61 0.79 0.78 0.78]';
% chi_max1  = [13.7 16.5 14.5 19.5 19.0 18.0]';
chi_max  = [13.7 16.7 26.8 28.2 24.8 16]'
zmax     = [6.2 10.5 6.1 22.0 21.0 22.0]';
lambda1  = [0.22 0.19 0.23 0.13 0.13 0.13]';
CSR = [0.2 0.25 0.3 0.5 0.4 0.3];
AlphaCI  = 0.05;

xg   = linspace(min(Dr)-0.01, max(Dr)+0.01, 500)';
R2fn = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);

%% ==================== FIGURA 1: AJUSTES + IC95% ====================
figure('Color','w','Name','Ajustes vs D_r (solo IC95%)');
set(gcf,'Units','pixels','Position',[80 120 1320 420]);

% ---------- (1) chi_max vs Dr  -> lineal ----------
subplot(1,3,1); hold on; box on; grid on;
mdl_chi = fitlm(Dr, chi_max);
[yhat_g, CI_g] = predict(mdl_chi, xg, 'Alpha',AlphaCI);
yhat_tr = predict(mdl_chi, Dr);
R2_fit  = R2fn(chi_max, yhat_tr);

fill([xg; flipud(xg)], [CI_g(:,1); flipud(CI_g(:,2))], [0.7 0.8 1.0], ...
    'EdgeColor','none', 'FaceAlpha',0.35);
scatter(Dr, chi_max, 42, 'filled');
plot(xg, yhat_g, 'b-', 'LineWidth',1.4);

xlabel('D_r (decimal)'); ylabel('\chi_{max}');
title(sprintf('\\chi_{max} vs D_r  |  R^2_{fit}=%.3f', R2_fit));
legend({'IC95% (media)','Datos','Ajuste (OLS)'}, 'Location','best');
xlim([min(xg) max(xg)]);
smartLabels(gca, Dr, chi_max, TCUI);

% ---------- (2) zmax vs Dr  -> no lineal c*x^k ----------
subplot(1,3,2); hold on; box on; grid on;
lm0   = fitlm(log(Dr), log(zmax));
b10   = exp(lm0.Coefficients.Estimate(1));
b20   = lm0.Coefficients.Estimate(2);
mod_z = @(b,x) b(1).*x.^b(2);

mdl_z      = fitnlm(Dr, zmax, mod_z, [b10 b20]);
[yhat_g,CI_g] = predict(mdl_z, xg, 'Alpha',AlphaCI);
yhat_tr  = predict(mdl_z, Dr);
R2_fit   = R2fn(zmax, yhat_tr);

fill([xg; flipud(xg)], [CI_g(:,1); flipud(CI_g(:,2))], [0.7 0.8 1.0], ...
    'EdgeColor','none', 'FaceAlpha',0.35);
scatter(Dr, zmax, 42, 'filled');
plot(xg, yhat_g, 'b-', 'LineWidth',1.4);

xlabel('D_r (decimal)'); ylabel('z_{max}');
title(sprintf('z_{max} vs D_r  |  R^2_{fit}=%.3f', R2_fit));
legend({'IC95% (media)','Datos','Ajuste (NLR)'}, 'Location','best');
xlim([min(xg) max(xg)]);
smartLabels(gca, Dr, zmax, TCUI);

% ---------- (3) lambda1 vs Dr  -> no lineal c*exp(-k x) ----------
subplot(1,3,3); hold on; box on; grid on;
mod_l1 = @(b,x) b(1).*exp(-b(2).*x);
b10 = max(lambda1); b20 = 1;
mdl_l1      = fitnlm(Dr, lambda1, mod_l1, [b10 b20]);
[yhat_g,CI_g] = predict(mdl_l1, xg, 'Alpha',AlphaCI);
yhat_tr  = predict(mdl_l1, Dr);
R2_fit   = R2fn(lambda1, yhat_tr);

fill([xg; flipud(xg)], [CI_g(:,1); flipud(CI_g(:,2))], [0.7 0.8 1.0], ...
    'EdgeColor','none', 'FaceAlpha',0.35);
scatter(Dr, lambda1, 42, 'filled');
plot(xg, yhat_g, 'b-', 'LineWidth',1.4);

xlabel('D_r (decimal)'); ylabel('\lambda_1');
title(sprintf('\\lambda_1 vs D_r  |  R^2_{fit}=%.3f', R2_fit));
legend({'IC95% (media)','Datos','Ajuste (NLR)'}, 'Location','best');
xlim([min(xg) max(xg)]);
smartLabels(gca, Dr, lambda1, TCUI);

set(findall(gcf,'Type','axes'),'FontName','Helvetica','FontSize',10);

%% ==================== FIGURA 2: DIAGNÓSTICOS ====================
figure('Color','w','Name','Diagnósticos: residuales studentizados');
set(gcf,'Units','pixels','Position',[80 80 1320 360]);

% chi_max
subplot(1,3,1); hold on; box on; grid on;
t = student_resid(mdl_chi);
plot_resid(gca, predict(mdl_chi,Dr), t, mdl_chi.Diagnostics.Leverage, TCUI);
xlabel('\hat{y} = \chi_{max}^{\,pred}'); ylabel('t_i'); title('\chi_{max}');

% zmax
subplot(1,3,2); hold on; box on; grid on;
t = student_resid(mdl_z);
plot_resid(gca, predict(mdl_z,Dr), t, mdl_z.Diagnostics.Leverage, TCUI);
xlabel('\hat{y} = z_{max}^{\,pred}');   ylabel('t_i'); title('z_{max}');

% lambda1
subplot(1,3,3); hold on; box on; grid on;
t = student_resid(mdl_l1);
plot_resid(gca, predict(mdl_l1,Dr), t, mdl_l1.Diagnostics.Leverage, TCUI);
xlabel('\hat{y} = \lambda_1^{\,pred}'); ylabel('t_i'); title('\lambda_1');

set(findall(gcf,'Type','axes'),'FontName','Helvetica','FontSize',10);
end  % ===== fin correlations_noPI() =====

% ===== Utilidades =====
function smartLabels(ax,x,y,labels)
    xl = xlim(ax); yl = ylim(ax);
    dx = 0.012*(xl(2)-xl(1));  dy = 0.020*(yl(2)-yl(1));
    pattern = [ 1  1; -1  1;  1 -1; -1 -1;  2 0; -2 0; 0  2; 0 -2];
    placed  = zeros(0,2);
    for k = 1:numel(x)
        pos = [x(k)+dx, y(k)+dy]; t = 0;
        while ~isempty(placed) && min(vecnorm(placed - pos,2,2)) < 1.2*max(dx,dy) && t < 24
            p   = pattern(mod(t,size(pattern,1))+1,:);
            pos = [x(k)+p(1)*dx, y(k)+p(2)*dy]; t = t+1;
        end
        text(ax, pos(1), pos(2), [' ' labels{k}], 'FontSize',8, ...
            'BackgroundColor','w', 'Margin',1, ...
            'VerticalAlignment','middle','Clipping','on');
        placed(end+1,:) = pos; %#ok<AGROW>
    end
end

function t = student_resid(mdl)
    r = mdl.Residuals.Raw;             
    h = mdl.Diagnostics.Leverage;      
    s = mdl.RMSE;                      
    t = r ./ (s .* sqrt(max(1e-12, 1 - h)));
end

function plot_resid(ax, yhat, t, lev, labels)
    scatter(yhat, t, 40, 'filled'); 
    yline(0,'k-'); yline(2,'--'); yline(-2,'--'); yline(3,':'); yline(-3,':');
    idxHi = abs(t) > 2 | lev > 2*mean(lev);
    arrayfun(@(k) text(yhat(k), t(k), [' ' labels{k}], ...
        'FontSize',8,'BackgroundColor','w','Margin',1), find(idxHi));
    ylim([-3.5 3.5]);
end


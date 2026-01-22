function J = fitnessfunction_karlsuhe_cyc2(x, y0, parms, nspb, path_info, p_exp, q_exp, eps1_exp)
% x = [cz, zmax, lambda1]
% Error medio adimensional en el plano eps1–q (cíclico).

% ---------- 0) Cotas físicas + salida rápida ----------
% lb = [ 300,  1.0, 0.10];
% ub = [2000, 40.0, 2.0];
lb = [500, 1.0, 0.10];  ub = [2000, 15, 1.0];

if any(x < lb) || any(x > ub)
    J = 1e6; return
end

% ---------- 1) Escalas y datos (columnas) ----------
q_exp    = q_exp(:);
eps1_exp = eps1_exp(:);
Ne       = numel(q_exp);
scale_q    = max(1, max(abs(q_exp)));
scale_eps1 = max(1, max(abs(eps1_exp)));

% ---------- 2) Parám. del modelo ----------
par = parms;           % copy-on-write
par(21) = x(1);        % cz
par(22) = x(2);        % zmax
par(27) = x(3);        % lambda1
par(33) = x(4);        % chimax

% (ajusta índices si tu layout difiere)

% ---------- 3) Correr driver ----------
try
    % Mantén 25,25 como pediste
    [SS, EE, ~, ~, ~] = updateModel_stress(y0, par, nspb, path_info, 25, 25);
catch
    J = 1e6; return
end

% ---------- 4) Serie del modelo (lazo eps1–q) ----------
q_mod    = SS(:,3) - SS(:,1);
eps1_mod = EE(:,3);
if isempty(q_mod) || any(~isfinite(q_mod)) || any(~isfinite(eps1_mod))
    J = 1e6; return
end

% ---------- 5) Re-muestreo del modelo a la malla experimental ----------
Nm = numel(q_mod);
if Nm ~= Ne
    s_m = linspace(0,1,Nm);
    s_e = linspace(0,1,Ne);
    q_mod    = interp1(s_m, q_mod,    s_e, 'pchip', 'extrap');
    eps1_mod = interp1(s_m, eps1_mod, s_e, 'pchip', 'extrap');
end

% **FORZAR COLUMNAS Y MISMA LONGITUD** (evita NxN por broadcasting)
q_mod    = reshape(q_mod   , [], 1);
eps1_mod = reshape(eps1_mod, [], 1);
q_exp    = reshape(q_exp   , [], 1);
eps1_exp = reshape(eps1_exp, [], 1);

% **VALIDAR** longitudes coherentes
Ne_m = numel(q_mod);
if Ne_m ~= Ne
    % como fallback, iguala al mínimo
    Nmin     = min(Ne, Ne_m);
    q_mod    = q_mod(1:Nmin);
    eps1_mod = eps1_mod(1:Nmin);
    q_exp    = q_exp(1:Nmin);
    eps1_exp = eps1_exp(1:Nmin);
    Ne       = Nmin;
end

% ---------- 6) Error medio adimensional (MSE) en eps1–q ----------
d_eps = (eps1_mod - eps1_exp) ./ scale_eps1;   % ahora son N×1
d_q   = (q_mod    - q_exp   ) ./ scale_q;
J = mean(d_eps.^2 + d_q.^2);

% ---------- 7) Plot en vivo (más liviano) ----------
persistent evalCount
if isempty(evalCount), evalCount = 0; end
evalCount = evalCount + 1;

if mod(evalCount, 10) == 1        % dibuja cada 10 evaluaciones
    % **Diezma solo para plot si N es grande**
    if Ne > 15000
        step = ceil(Ne/15000);    % ~<= 15k puntos en pantalla
        idx  = 1:step:Ne;
    else
        idx = 1:Ne;
    end

    figure(1);
    plot(eps1_exp(idx), q_exp(idx), 'k-', 'LineWidth', 0.5); hold on
    plot(eps1_mod(idx), q_mod(idx), 'r-', 'LineWidth', 0.5); hold off
    xline(0,'--','color',[.5 .5 .5],'LineWidth',0.5);
    yline(0,'--','color',[.5 .5 .5],'LineWidth',0.5);
    xlabel('\epsilon_1'); ylabel('q [kPa]'); grid on; box on
    title(sprintf('J=%.3e | cz=%.0f, z_{max}=%.2f, \\lambda_1=%.3f (N=%d)', ...
          J, x(1), x(2), x(3), Ne));
%     pause(0.03);
end
end
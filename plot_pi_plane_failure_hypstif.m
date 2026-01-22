function out = plot_pi_plane_failure_hypstif(opts)
% plot_pi_plane_failure_hypstif
% π-plane locus centered & "flat" using the same a1->a0->F->ag framework as hyp_stif.txt.
% Produces ONE figure with LaTeX labels and principal axes arrows extending beyond the surface.

if nargin < 1, opts = struct(); end
opts = setDefaults(opts);

% -----------------------------
% Parameters (Toyoura default)
% -----------------------------
phi_c  = opts.phi_c_deg * pi/180;
eta    = opts.eta;                 % generalized failure parameter (0 => circle)
nPts   = opts.nPts;

% a_DP (same constant used in hyp_stif for ag without F)
aDP = sqrt(3)/2 * (3 - sin(phi_c)) / (sqrt(2)*sin(phi_c));

% π-plane angular parametrization (diagonal deviatoric directions)
alpha = linspace(0, 2*pi, nPts).';
ca = cos(alpha);  sa = sin(alpha);

% Deviatoric orthonormal basis (diagonal states)
% E1 = diag([ 2, -1, -1 ]) / sqrt(6)
% E2 = diag([ 0,  1, -1 ]) / sqrt(2)
n11 = (2/sqrt(6))*ca;
n22 = (-1/sqrt(6))*ca + (1/sqrt(2))*sa;
n33 = (-1/sqrt(6))*ca - (1/sqrt(2))*sa;

% ---------------------------------------------------------
% Solve the implicit critical locus:  r = ag(r, alpha)
% where ag = aDP * F,  F = exp(-eta * a0),  a0 = ln(a1)
% and a1 is computed from invariants (Yao-type) as in hyp_stif.
% ---------------------------------------------------------
r_gen = solve_r_fixed_point(n11, n22, n33, aDP, eta, opts);

% Isotropic reference (eta = 0 => F = 1 => r = aDP)
r_iso = aDP * ones(size(alpha));

% Rotate the plot to match the paper-like orientation (sigma11 up)
rot = opts.rotation; % default pi/2
x_gen = r_gen .* cos(alpha + rot);
y_gen = r_gen .* sin(alpha + rot);
x_iso = r_iso .* cos(alpha + rot);
y_iso = r_iso .* sin(alpha + rot);

% -----------------------------
% Plot (single centered figure)
% -----------------------------
figure('Name','\pi-plane failure locus (centered)','Color','w');
ax = axes(); hold(ax,'on'); box(ax,'on'); axis(ax,'equal');

% Curves
plot(ax, x_gen, y_gen, '-',  'LineWidth', 2.2);
if opts.showIsotropic
    plot(ax, x_iso, y_iso, '--', 'LineWidth', 1.6);
end

% Symmetric limits (this is what visually "centers" the plot)
rSurf = max([max(r_gen), max(r_iso)]);
rAxis = opts.axisRayFactor * rSurf;   % arrows go beyond surface
rLim  = opts.axisLimFactor * rSurf;   % view window

xlim(ax, [-rLim, rLim]);
ylim(ax, [-rLim, rLim]);

% Principal axes arrows (σ11, σ22, σ33) "coming out" of the surface
angAxes = rot + [0, 2*pi/3, 4*pi/3]; % σ11 along E1, then +120°, +240°
labs = {'$\hat{\sigma}_{11}$', '$\hat{\sigma}_{22}$', '$\hat{\sigma}_{33}$'};
for i = 1:3
    quiver(ax, 0, 0, rAxis*cos(angAxes(i)), rAxis*sin(angAxes(i)), 0, ...
        'k', 'LineWidth', 1.1, 'MaxHeadSize', 0.08);
    text(ax, 1.06*rAxis*cos(angAxes(i)), 1.06*rAxis*sin(angAxes(i)), labs{i}, ...
        'Interpreter','latex', 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
        'FontSize', opts.fontSize);
end

% Origin marker
plot(ax, 0, 0, 'k+', 'LineWidth', 1.2, 'MarkerSize', 8);

% Aesthetics (LaTeX ticks/labels)
ax.TickLabelInterpreter = 'latex';
ax.FontSize  = opts.fontSize;
ax.LineWidth = 1.0;
grid(ax, 'on');

xlabel(ax, '$x\;(\hat{\boldsymbol{\sigma}}^{*}:\mathbf{E}_1)$', 'Interpreter','latex', 'FontSize', opts.fontSize+1);
ylabel(ax, '$y\;(\hat{\boldsymbol{\sigma}}^{*}:\mathbf{E}_2)$', 'Interpreter','latex', 'FontSize', opts.fontSize+1);

ttl = sprintf('\\pi-plane failure locus (hyp\\_stif framework), \\;\\eta = %.3g', eta);
title(ax, ttl, 'Interpreter','latex', 'FontSize', opts.fontSize+2);

if opts.showIsotropic
    legend(ax, {'Generalized (hyp\_stif)','Isotropic ($\eta=0$)'}, ...
        'Interpreter','latex', 'Location','northeast');
else
    legend(ax, {'Generalized (hyp\_stif)'}, 'Interpreter','latex', 'Location','northeast');
end

% Output
out = struct();
out.alpha = alpha;
out.r_gen = r_gen;
out.r_iso = r_iso;
out.xy_gen = [x_gen, y_gen];
out.xy_iso = [x_iso, y_iso];
out.aDP   = aDP;
out.eta   = eta;
out.phi_c_deg = opts.phi_c_deg;

end

% ==========================================================
% Defaults
% ==========================================================
function opts = setDefaults(opts)
def = struct();
def.phi_c_deg     = 34;      % Toyoura-ish
def.eta           = 0.1;    % try 0.2–0.8 (0 => circle)
def.nPts          = 721;

% Fixed-point solver controls (vectorized)
def.maxIter       = 60;
def.relax         = 0.75;    % under-relaxation for robustness
def.tol           = 5e-12;
def.rMin          = 1e-6;
def.rMax          = 1.50;    % safety clamp (dimensionless)
def.sqrtEps       = 1e-14;
def.a1Min         = 1e-10;

% Plot controls
def.showIsotropic = true;
def.rotation      = pi/2;    % rotate so σ11 is upward
def.axisRayFactor = 1.18;    % arrows beyond surface
def.axisLimFactor = 1.30;    % symmetric window beyond arrows
def.fontSize      = 12;

fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opts, fn{k}) || isempty(opts.(fn{k}))
        opts.(fn{k}) = def.(fn{k});
    end
end
end

% ==========================================================
% Vectorized fixed-point solver: r = aDP * a1(r,alpha)^(-eta)
% ==========================================================
function r = solve_r_fixed_point(n11, n22, n33, aDP, eta, opts)

% Start from isotropic radius
r = aDP * ones(size(n11));

if eta == 0
    return; % exactly circular
end

for it = 1:opts.maxIter
    r_old = r;

    % Diagonal normalized stress: T = I/3 + r*n
    T11 = 1/3 + r .* n11;
    T22 = 1/3 + r .* n22;
    T33 = 1/3 + r .* n33;

    % Invariants (DIAGONAL ONLY, consistent with your π-plane construction)
    I1 = T11 + T22 + T33;  % should be 1
    I2 = 0.5 * ( (T11.^2 + T22.^2 + T33.^2) - I1.^2 );
    I3 = T11 .* T22 .* T33;

    % Yao-type a1 (same algebraic structure as in hyp_stif)
    a11 = 2 .* I1;

    num = (I1.*I2 - I3) .* (I1.^2 - 3.*I2);
    den = (I1.*I2 - 9.*I3);

    arg12 = num ./ den;
    arg12 = max(arg12, opts.sqrtEps);

    arg13 = I1.^2 - 3.*I2;
    arg13 = max(arg13, opts.sqrtEps);

    a12 = 3 .* sqrt(arg12);
    a13 = sqrt(arg13);

    denom = (a12 - a13);
    % Guard against division by ~0 and sign flips
    denom(abs(denom) < 1e-14) = sign(denom(abs(denom) < 1e-14)) * 1e-14;

    a1 = a11 ./ denom;
    a1 = max(a1, opts.a1Min);

    % F = exp(-eta * ln(a1)) = a1^(-eta)
    F  = a1 .^ (-eta);
    ag = aDP .* F;

    % Fixed-point update with relaxation
    r_new = ag;
    r = opts.relax * r_new + (1 - opts.relax) * r_old;

    % Clamp
    r = max(r, opts.rMin);
    r = min(r, opts.rMax);

    % Convergence
    if max(abs(r - r_old)) < opts.tol
        break;
    end
end

end

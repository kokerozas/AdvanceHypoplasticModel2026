function fig5_flat_toyoura_pretty()
% Fig. 5(A)-style "flat + centered" π-plane plot (Toyoura), improved styling:
% - LaTeX labels
% - axes extend beyond the locus (visual)
% - clean, paper-like layout (single figure)

% -------------------------
% Toyoura parameters (Table 1)
% -------------------------
phi_c_deg = 34;      % [deg]
Fin       = 0.45;    % initial fabric magnitude
lambda    = 0.36;    % anisotropic strength scaling
N         = 1600;    % angular resolution

% -------------------------
% c1,c2 from phi_c (Eq.10a-10b)
% -------------------------
phi_c = deg2rad(phi_c_deg);
c1 = sqrt(3/8) * (3 - sin(phi_c)) / sin(phi_c);
c2 = (3/8)     * (3 + sin(phi_c)) / sin(phi_c);

% -------------------------
% Unit deviatoric directions in principal space (diagonal, tr(n)=0, ||n||=1)
% n(alpha) = cos(alpha)E1 + sin(alpha)E2 -> components (n11,n22,n33)
% -------------------------
alpha = linspace(0,2*pi,N).';
ca = cos(alpha); sa = sin(alpha);

n11 = (2/sqrt(6))*ca;
n22 = (-1/sqrt(6))*ca + (1/sqrt(2))*sa;
n33 = (-1/sqrt(6))*ca - (1/sqrt(2))*sa;

tr2 = n11.^2 + n22.^2 + n33.^2;    % = 1 (numerical)
tr3 = n11.^3 + n22.^3 + n33.^3;

cos3 = -sqrt(6) .* tr3 ./ (tr2.^(3/2));
g    = abs(1 + cos3);

% -------------------------
% Initial fabric direction (cross-anisotropic): nF = diag(2,-1,-1)/sqrt(6)
% F_in = Fin * nF
% -------------------------
nF11 =  2/sqrt(6);
nF22 = -1/sqrt(6);
nF33 = -1/sqrt(6);

A  = Fin .* (nF11*n11 + nF22*n22 + nF33*n33); % A(alpha)=F:n

% -------------------------
% a1(alpha) from: c2*g*a1^2 + c1*a1 - rhs = 0 (positive root)
% isotropic rhs=1, anisotropic rhs=exp(lambda*(A-1))
% -------------------------
rhs_iso  = ones(size(alpha));
rhs_anis = builtin('exp', lambda .* (A - 1));

a1_iso  = solve_a1(c1,c2,g,rhs_iso);
a1_anis = solve_a1(c1,c2,g,rhs_anis);

% Principal deviatoric components of \hat{s}^*:
d1_iso = a1_iso .* n11;  d2_iso = a1_iso .* n22;  d3_iso = a1_iso .* n33;
d1_an  = a1_anis.* n11;  d2_an  = a1_anis.* n22;  d3_an  = a1_anis.* n33;

% -------------------------
% 2D Fig.5-style tri-axial 120° embedding
% axes: \hat{s}^*_{11} (up), \hat{s}^*_{22} (down-left), \hat{s}^*_{33} (down-right)
% -------------------------
b11 = [ 0;  1];
b22 = [-sqrt(3)/2; -1/2];
b33 = [ sqrt(3)/2; -1/2];

Piso = (d1_iso.' .* b11) + (d2_iso.' .* b22) + (d3_iso.' .* b33); % 2xN
Pan  = (d1_an.'  .* b11) + (d2_an.'  .* b22) + (d3_an.'  .* b33); % 2xN
Piso = [Piso Piso(:,1)];
Pan  = [Pan  Pan(:,1)];

% -------------------------
% Figure styling controls
% -------------------------
axisScale      = 1.35; % >1 makes axes "exit" the surface (longer rays)
frameMargin    = 1.10; % additional margin beyond axis length
labelOffset    = 1.06; % label distance beyond ray tip
labelBg        = 'w';  % small white patch behind labels
labelMargin    = 2;    % px

% ray length based on locus size
rLocus = max(abs([Pan(:); Piso(:)]));
rAxis  = axisScale * rLocus;
rFrame = frameMargin * rAxis;

% -------------------------
% Plot (ONE FIGURE)
% -------------------------
figure('Name','Fig.5(A)-style (pretty)'); clf;
ax = axes(); hold(ax,'on'); box(ax,'off'); axis(ax,'equal');
set(ax,'XTick',[],'YTick',[],'Visible','off'); % paper-like: no ticks, no frame

% curves
plot(Pan(1,:),  Pan(2,:),  '-',  'LineWidth', 2.2);
plot(Piso(1,:), Piso(2,:), '--', 'LineWidth', 1.6);

% axes rays (extend beyond surface)
drawRay(ax, b11, rAxis);
drawRay(ax, b22, rAxis);
drawRay(ax, b33, rAxis);

% centered frame
xlim(ax, [-rFrame rFrame]);
ylim(ax, [-rFrame rFrame]);

% LaTeX axis labels (placed beyond ray tips)
placeLabel(ax, b11, rAxis*labelOffset, '$\hat{\sigma}_{11}^\*$' , labelBg, labelMargin);
placeLabel(ax, b22, rAxis*labelOffset, '$\hat{\sigma}_{22}^\*$' , labelBg, labelMargin);
placeLabel(ax, b33, rAxis*labelOffset, '$\hat{\sigma}_{33}^\*$' , labelBg, labelMargin);

% Legend (also LaTeX)
lgd = legend(ax, {'Anisotropic criterion','Isotropic criterion'}, ...
    'Location','northeast');
set(lgd,'Interpreter','latex','Box','off');

% Title (LaTeX)
ttl = title(ax, 'Fig.~5(A)-style: $\pi$-plane locus (Toyoura)', 'Interpreter','latex');
set(ttl,'FontWeight','normal');

end

% =====================================================================
% Helpers
% =====================================================================

function a1 = solve_a1(c1,c2,g,rhs)
A = c2 .* g;  B = c1;
a1 = zeros(size(rhs));
lin = abs(A) < 1e-14;
qua = ~lin;

if any(lin)
    a1(lin) = rhs(lin) ./ B;
end
if any(qua)
    Aq = A(qua); rq = rhs(qua);
    disc = B.^2 + 4.*Aq.*rq;
    a1(qua) = (-B + sqrt(max(disc,0))) ./ (2.*Aq);
end
end

function drawRay(ax, b, L)
% Draws an axis ray from origin to L*b, with arrowhead-ish look.
p0 = [0;0];
p1 = L*b(:);
plot(ax, [p0(1) p1(1)], [p0(2) p1(2)], 'k-', 'LineWidth', 1.2);

% Add a small arrowhead by two short segments (robust, no quiver dependencies)
ah = 0.06*L;    % arrowhead size
ang = atan2(b(2), b(1));
phi = deg2rad(22);

p2 = p1 - ah*[cos(ang-phi); sin(ang-phi)];
p3 = p1 - ah*[cos(ang+phi); sin(ang+phi)];
plot(ax, [p1(1) p2(1)], [p1(2) p2(2)], 'k-', 'LineWidth', 1.2);
plot(ax, [p1(1) p3(1)], [p1(2) p3(2)], 'k-', 'LineWidth', 1.2);
end

function placeLabel(ax, b, L, txt, bg, marginPx)
% Places LaTeX label at L*b, aligned with ray direction.
p = L*b(:);
ha = 'center'; va = 'middle';

% slight alignment tweak by quadrant (readability)
ang = atan2(b(2), b(1));
if ang > 0.25*pi && ang < 0.75*pi
    va = 'bottom';
elseif ang < -0.25*pi && ang > -0.75*pi
    va = 'top';
end
if abs(ang) < 0.25*pi
    ha = 'left';
elseif abs(ang) > 0.75*pi
    ha = 'right';
end

t = text(ax, p(1), p(2), txt, ...
    'Interpreter','latex', ...
    'HorizontalAlignment',ha, ...
    'VerticalAlignment',va, ...
    'FontSize',12);

% background patch for clean overlap handling
set(t, 'BackgroundColor', bg, 'Margin', marginPx);
end


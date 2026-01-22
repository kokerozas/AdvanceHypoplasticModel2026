% plotting section
%
% sand hypoplastic constitutive model

array_sizes = size(INV_E);
nrow = array_sizes(1);

% p-e in oed compression

xx=INV_S(1:nrow,1);
% xx = (SS(:,3)+ SS(:,1))/2;
% xx=SS(1:nrow,3);
yy=HARD(1:nrow,1);

% subplot(1,2,1)
semilogx(xx,yy,'r-');hold on
xlabel(' $\sigma_1$ (kPa) ','interpreter','latex')
ylabel(' Void ratio, e (-)','interpreter','latex')
legend('Model','Experiments');
ylim([0.65 1.05])
xlim([1.0 1e3])
% Personalización de ejes para un estilo minimalista
% set(gca, 'FontSize', 10, 'FontName', 'Times New Roman', 'LineWidth', 0.8);  % Establecer tamaño de fuente y grosor de las líneas
% set(gca, 'Box', 'on', 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2]);  % Colores suaves para los ejes
% legend(h2, 'Model', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 10, 'Box', 'off');

% Activar cuadrícula sutil
% grid off;
% set(gca, 'LineWidth', 1.1);  % Establecer el grosor del borde del gráfico a 1.0

% Etiquetas y título con un estilo sutil en LaTeX
% xlabel(' p (kPa)', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
% ylabel('e (-)', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
% title('Pressure vs Void Ratio', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'normal');


% grid on





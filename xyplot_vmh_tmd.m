    % plotting section
%
% sand hypoplastic constitutive model

array_sizes = size(INV_E);
nrow = array_sizes(1);
%% TMU
% p-q
% subplot(2,2,3)
% % xx = (SS(:,3)+SS(:,1))/2;
% yy = SS(:,3)-SS(:,1);
% xx= INV_S(1:nrow,1);
% % yy=INV_S(1:nrow,2);
% plot(xx,yy,'r-','linewidth',0.5); hold on
% xlabel('$p$ (kPa)', 'interpreter', 'latex')
% ylabel('$q$ (kPa)','interpreter','latex')
% % xticks(0:100:400);    xlim([0,400]);
% % yticks(-300:150:450); ylim([-300,450]);
% xticks(0:100:600);    xlim([0,600]);
% yticks(-400:200:800); ylim([-400,800]);
% xline(0, '--','color',[0.5 0.5 0.5], 'LineWidth', 0.5); % Línea horizontal en q = 50 kPa
% yline(0, '--','color',[0.5 0.5 0.5], 'LineWidth', 0.5); % Línea horizontal en q = 50 kPa
% text(-0.18, 1.05, '(C)', ...
%     'Units', 'normalized', ...
%     'FontWeight', 'normal', ...
%     'FontSize', 11, ...
%     'Interpreter', 'Latex', ...
%     'HorizontalAlignment', 'left', ...
%     'VerticalAlignment', 'top', ...
%     'Clipping', 'off')
% 
% 
% % eps1-q
% subplot(2,2,4)
% % xx=INV_E(1:nrow,2);
% % yy = SS(:,3)-SS(:,1);
% xx = EE(:,3);
% yy= INV_S(1:nrow,2);
% 
% plot(xx,-yy,'r-','linewidth',0.5);hold on
% ylabel('$q$ (kPa)', 'interpreter', 'latex')
% xlabel('$\varepsilon_1$ (-)', 'interpreter', 'latex')
% xticks(-0.1:0.05:0.1);xlim([-0.1,0.1]);
% yticks(-500:500:1000); ylim([-500,1000]);
% xline(0, '--','color',[0.5 0.5 0.5], 'LineWidth', 0.5); % Línea horizontal en q = 50 kPa
% yline(0, '--','color',[0.5 0.5 0.5], 'LineWidth', 0.5); % Línea horizontal en q = 50 kPa
% text(-0.18, 1.05, '(D)', ...
%     'Units', 'normalized', ...
%     'FontWeight', 'normal', ...
%     'FontSize', 11, ...
%     'Interpreter', 'Latex', ...
%     'HorizontalAlignment', 'left', ...
%     'VerticalAlignment', 'top', ...
%     'Clipping', 'off')

% set(gcf, 'PaperUnits', 'centimeters')
% set(gcf, 'PaperSize', [37 22])
% set(gcf, 'PaperPosition', [0 0 37 22])

%% TMD
% % % epss-q plot
subplot(1,2,1)
xx= EE(:,3) ;
yy= (INV_S(:,2)) ;
% xline(0.0557,':');
% plot(xx,yy,'color',[0.72,0.27,1.00]); hold on
plot(xx,yy,'k-');hold on 
xlabel('$\varepsilon_1$ (-)', 'interpreter', 'latex')
ylabel('$ q$ (kPa)','interpreter','latex')
xlim([0 0.3]);xticks(0:0.1:0.3); 
% ylim([0 1600]);yticks(0:400:1600); 

subplot(1,2,2)
xx= INV_E(:,2) ;
yy= INV_E(:,1) ;
% xline(0.0557,':');
% plot(xx,yy,'color',[0.72,0.27,1.00]); hold on
plot(xx,yy,'k-');hold on 
xlabel('$\varepsilon_1$ (-)', 'interpreter', 'latex')
ylabel('$ q$ (kPa)','interpreter','latex')


xlim([0 0.3]);xticks(0:0.1:0.3); 
% ylim([0 1600]);yticks(0:400:1600); 



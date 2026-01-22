% plotting section
%
% sand hypoplastic constitutive model

array_sizes = size(INV_E);
nrow = array_sizes(1);


%clf
% 
% epss-q plot
subplot(2,2,1)
xx=INV_E(1:nrow,2);
yy=INV_S(1:nrow,2);

plot(xx,yy,'r-')
xlabel('deviatoric strain')
ylabel('deviator stress q [kPa]')
grid on
hold on

% p-q plot
subplot(2,2,2)
xx=INV_S(1:nrow,1);
yy=INV_S(1:nrow,2);
plot(xx,yy,'r-')
xlabel('mean effective stress p [kPa]')
ylabel('deviator stress q [kPa]')
grid on
hold on


% epsv-epss plot
subplot(2,2,3)
xx=INV_E(1:nrow,2);
yy=INV_E(1:nrow,1);
plot(xx,yy,'r-')

xlabel('deviatoric strain')
ylabel('volumetric strain')
grid on
hold on



% eps_v-p plot

subplot(2,2,4)
xx=INV_E(1:nrow,1);
yy=INV_S(1:nrow,2);

plot(xx,yy,'r-')
xlabel('volumetric strain')
ylabel('q [kPa]')
grid on
hold on

% figure(2)
%clf

% c-eps_s plot
% 
% xx=INV_E(1:nrow,2);
% yy=HARD(1:nrow,1);
% 
% plot(xx,yy,'r-')
% ylabel('void ratio [-]')
% xlabel('deviatoric strain')
% grid on
% hold on


%% 

% 
% figure(2)
% subplot(1,2,2)
% xx=INV_E(1:nrow,2);
% yy=INV_E(1:nrow,1);
% plot(xx,yy,'b-','Linewidth',1.1)
% 
% ylabel('Deformación Volumetrica', 'Fontsize' , 12 , 'Fontname' , 'Calisto MT')
% xlabel('Deformación Axial', 'Fontsize' , 12 , 'Fontname' , 'Calisto MT')
% grid on
% hold on
% 
% set(gca,'Fontname','Calisto MT')
% grid on
% 
% 
% subplot(1,2,1)
% xx=INV_E(1:nrow,2);
% yy=INV_S(1:nrow,2);
% 
% plot(xx,yy,'b-','Linewidth',1.1)
% xlabel('Deformación Axial', 'Fontsize' , 12 , 'Fontname' , 'Calisto MT')
% ylabel('Desviador de Tensiones q (kPa)', 'Fontsize' , 12 , 'Fontname' , 'Calisto MT')
% grid on
% hold on
% 
% set(gca,'Fontname','Calisto MT')
% grid on

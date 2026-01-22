function [Mep,Hdel,d_strain,d_strain2] = M_istr(Lep,Nep,deps,istrain,parms,Strain_total,FF,Strain_total2)

% computes matrix M of the intergranular strain concept by Niemunis and Herle (1997)

M2=[1,  0,  0,  0,    0,    0
   0,  1,  0,  0,    0,    0
   0,  0,  1,  0,    0,    0 
   0,  0,  0,  0.5,  0,    0
   0,  0,  0,  0,    0.5,  0
   0,  0,  0,  0,    0,    0.5];

I=[1,  0,  0,  0,  0,  0
   0,  1,  0,  0,  0,  0
   0,  0,  1,  0,  0,  0
   0,  0,  0,  1,  0,  0
   0,  0,  0,  0,  1,  0
   0,  0,  0,  0,  0,  1];

mR= parms(12);
mT= parms(13);
Rparam= parms(14) ;
betar= parms(15);
chimax=parms(33);

% chi  = parms(16) + 2*parms(16)*sin(Strain_total*pi*0.8333); % depende de los ciclos
% betar= parms(15) - FF*0.3*parms(15)*Strain_total; % depende de los ciclos
% chi= (parms(16) + 1.55*parms(16)*max(Strain_total,0)); % depende de los ciclos

% if p > FF*2
% betar= parms(15) - 0.7*parms(15)*(exp(-5*exp(-(Strain_total)*50))  - exp(-5*exp(-(Strain_total-0.8)*30)) );

% if Strain_total == 1
%     betar = parms(15);
% else
% betar= parms(15) - 0.5*parms(15)*Strain_total; 
% end

%     Rparam  = (parms(14) + 0.5*parms(14)*(exp(-5*exp(-(Strain_total-0.1)*15)))); % depende de los ciclos
% chi  = (parms(16) + 1.55*parms(16)*Strain_total*Factor2); % depende de los ciclos
% chi    = (parms(16) + (2.0*parms(16))*max(Strain_total,0)); % depende de los ciclos
% chi    = parms(16) + max(Strain_total,0)*(16.7 - parms(16));
% chi    = parms(16) + max(Strain_total,0)*(1.0 - parms(16));


% chi    = parms(16) + max(Strain_total,0)*(15 - parms(16));

% chi    = parms(16);
%Poblete 2016
chi    = (parms(16) + (Strain_total *(chimax - parms(16))));
% chi    = (parms(16) + (Strain_total *(12 - parms(16))));
%%.
% Al inicio de la carga cíclica (fabric intacto, $Zdd$ bajo): $FF \approx 1$, entonces chi aumenta según la fórmula casi linealmente con Strain_total. La estructura todavía puede degradarse más, por lo que se permite que chi suba ...
...Tras muchos ciclos: $Zdd$ se acerca a $z_{max}$ (se ha agotado la capacidad de cambio estructural), entonces $FF \to 0$. En ese límite, la fórmula de $\chi$ vuelve prácticamente a $\chi = \text{chi}_0$ (pues $FF \cdot$Strain_total se anula)...
... Esto significa que ya no se incrementa más chi, es decir, la degradación de rigidez por fabric se ha saturado. El modelo interpreta que la estructura se ha remodelado completamente y no hay un efecto adicional por más ciclos.

% chinew = p00 + p10*x + p01*y
% p00 =      -5.421  (-8.632, -2.209)
% p10 =       8.231  (4.175, 12.29)
% p01 =      0.0447  (0.03274, 0.05666)

istrain_norm=sqrt(istrain'*M2*istrain);
TINY=1.e-20;
istrain_dir=zeros(6,1);

if  (istrain_norm>TINY)
     istrain_dir=istrain/istrain_norm;
end

istraindir_D=istrain_dir'*M2*deps;
rho=istrain_norm/Rparam;

if rho > 1
    rho = 1;
end

LijmnSmn=Lep*istrain_dir;
LijmnSmnSkl=LijmnSmn*istrain_dir';
NijSkl=Nep*istrain_dir'; % dyadic product properties

Mep=(mT*rho^chi+mR*(1-rho^chi))*Lep;
% Mep=mT*Lep; %monotonic condition

if(istraindir_D>0) 
  Mep=Mep+(1-mT)*rho^chi*LijmnSmnSkl-rho^chi*NijSkl;
  Hdel=I-rho^betar*M2*istrain_dir*istrain_dir';
else
  Mep=Mep+(mR-mT)*rho^chi*LijmnSmnSkl;
  Hdel=I;
end



dev_strain_rate_norm = sqrt(deps'*I*deps);
yh = (rho^chi*Nep)'*deps;
%deps := incremento de deformación aplicado
% d_strain = 5.0*(1 - yh - Strain_total) * dev_strain_rate_norm;
d_strain = 20*(1 -(rho^chi*Nep)'*deps - Strain_total) *dev_strain_rate_norm;

d_strain2 = 15 *(1 - Strain_total2)*dev_strain_rate_norm;
% d_strain2 = 15.3 *(1 - Strain_total2)*dev_strain_rate_norm;
% d_strain2 = 1.0 *(1 - Strain_total2)*dev_strain_rate_norm;

% if Strain_total < 0
%     d_strain = 0;
% end




global log_internal_vars log_index
log_internal_vars.d_strain(log_index) = d_strain;
log_internal_vars.Strain_total(log_index) = Strain_total;
log_internal_vars.Strain_total2(log_index) = Strain_total2;
log_internal_vars.d_strain2(log_index) = d_strain2;
log_internal_vars.chi(log_index) = chi;
log_internal_vars.yh(log_index) = yh;
% disp(['Chi actual: ', num2str(chi), ' | FF: ', num2str(FF), ' | Strain_total: ', num2str(Strain_total), ' | d_strain: ', num2str(d_strain), ' | d_strain2: ', num2str(d_strain2)]);

end

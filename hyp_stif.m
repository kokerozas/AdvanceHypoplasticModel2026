function [Lep,Nep,Hep,z_dot, l_dot, Se, dE,p,FF,dZZ,d_strain2,l_dot2] = hyp_stif(y,parms,deps)
% 
% HYP_STIF: computes hypoplastic tensors L and N
%

M=[1,  0,  0,  0,  0,  0
   0,  1,  0,  0,  0,  0
   0,  0,  1,  0,  0,  0
   0,  0,  0,  2,  0,  0
   0,  0,  0,  0,  2,  0
   0,  0,  0,  0,  0,  2];

I=[1,  0,  0,  0,     0,     0
   0,  1,  0,  0,     0,     0
   0,  0,  1,  0,     0,     0 
   0,  0,  0,  0.5,   0,     0
   0,  0,  0,  0,     0.5,   0
   0,  0,  0,  0,     0,   0.5];

kron_delta=[1,  1,  1,  0,  0,  0];

ny = max(size(y));
sig = y(7:12,1);
DD  = y(1:6,1);
qint = y(13:ny,1);
% Strain_total = qint(2:7,1);


%Ftensor
zint = y(20:25,1);
%Semifluized state
lint = y(26,1);
lint2= y(31,1);

% En   = y(27,1);
% Strain_cum = y(28,1);

%Magnitud del cambio del tensor zint
Zdd = y(29,1);

% add cohesion
sig = sig + kron_delta'*parms(18)/tand(parms(4));

e=qint(1);
[p,qdev,z]   = inv_s(sig);
[dev,des,de] = inv_e(deps);
[DDv,DDs,~]  = inv_e(DD);
[Dv,Ds,~]    = inv_e(deps);


T=sig;
T_str=T/(3*p); % normal stress normalized 
% Definición de tensores 
T_star = sig-p*kron_delta';         %parte desviatórica
T_str_star = T_str-kron_delta'/3;   %desviatorica corregida

[p_str,qdev_str,z_str] = inv_s(T_str);

TT_str=stress_mult(T_str,T_str);
TT=stress_mult(T_str_star,T_str_star);
TTT=stress_mult(TT,T_str_star);

Tstarnorm = sqrt(T_star'*M*T_star); %norma tensor desviatorico
Tstrstarnorm = sqrt(T_str_star'*M*T_str_star); %norma tensión desviatorica corregida

phic        = parms(4)/180*3.14159265358979;
hs       	= parms(5);
nparam     	= parms(6);
ed0 		= parms(7);
ec0	     	= parms(8);
ei0	   	    = parms(9);
alpha	   	= parms(10);
beta	   	= parms(11);
hs_reverse  = parms(17);

tanpsi=sqrt(3.)*Tstrstarnorm;
TTtrace=stress_trace(TT);
if stress_trace(TT)==0
  TTtrace=0.0000000000000001;
end
%% Barotropy, Picknotropy and Volumetric Factors

cos3theta=-sqrt(6.)*stress_trace(TTT)/(TTtrace^1.5);
%Matsuoka Nakai critical state surface criterium
% F1=sqrt((tanpsi*tanpsi/8)+((2-tanpsi*tanpsi)/(2+sqrt(2.)*tanpsi*cos3theta)))-(tanpsi/(2*sqrt(2.)));

% critic friction angle 
% phic = phic;

% phic = phic*(1-beta) + beta*phic*y(30,1);
phic = phic*( (1-beta) + beta*y(30,1));


%Void ratios definition
expe=-1*(3*p/hs)^nparam;
ei=ei0*exp(expe);
ec=ec0*exp(expe);
ed=ed0*exp(expe);

% smalldiff=0.0001;
% re=(e-ed)/(ec-ed);
% fd=re^alpha;

% if e<=(1+smalldiff)*ed
%   e = ed;
%   disp('e lower than ed');
% end
% 
% if e>=(1+smalldiff)*ei
%   e = ei;
%   disp('e lower than ed');
% end

r = sqrt(T_star'*M*T_star)/p;
% rd   = (ei-parms(32))/(ei-ed);

% I_pe = (ei/e)*exp((ec-e))/((1+r));
% I_pe = 1.05*rd^(rd*alpha)*(ei/e)*exp((ec-e)^3)/((1+r));
% I_pe = rd^(rd*(1-alpha))*(ei/e)^(1+rd)*exp((ec-e))/((1+r)^beta);
% I_pe = 1;
% new expression
I_pe = ((ei/e)*exp((ec-e))/((1+r)))^exp(ec-e); 



% %stress invariants Yao et al., 2020
% if qdev > 1e-6
I1 = T_str(1) + T_str(2) + T_str(3);
I2 = 0.5*(T_str'*M*T_str - I1^2);
I3 = T_str(1)*T_str(2)*T_str(3);

a11= 2*I1;
a12= 3*sqrt((I1*I2-I3)*(I1^2-3*I2)/(I1*I2-9*I3));
a13= sqrt(I1^2-3*I2); % q expressed in form of stress invariants
a1 = a11/(a12-a13); %
a0 = log(a1);
% a01= max(real(1/exp(a0*0.2)),1.0); %valor móvil (0.0,0.5 and 1.0) 
% a01= real(1/exp(a0*0.2)); %valor móvil (0.0,0.5 and 1.0) 
a01= real(1/exp(a0*0.0)); %valor móvil (0.0,0.5 and 1.0) 



% r_parms = parms(19)*a01; %Material constant Ki/Gi * 
r_parms = parms(19); %Material constant Ki/Gi


exp1 = nparam*e*(3*p/hs)^nparam;
% exp1 =(3*p/hs)^nparam;
% exp2 = nparam*e*(3*p/(parms(17)*hs))^(nparam);

ag = sqrt(3)/2*(3-sin(phic))/(sqrt(2)*sin(phic));
ag=ag*a01; %strength parameters *a01

% ag = sqrt(3)/2*(3-sin(phic))/(sqrt(2)*sin(phic)); 
% f_a= 1+1/3*ag^2-sqrt(3)/3; %Herle and Kolymbas, 2004
f_a= 1+1/3*ag^2-sqrt(3)/3*ag; %Herle and Kolymbas, 2004
f_v= 3/2*r_parms - f_a;
fd   = (e/ec)^(alpha);

% %Original parameters
% a=sqrt(3.)*(3-sin(phic))/(2*sqrt(2.)*sin(phic));
% fs=(hs/nparam)*(ei/e)^beta*((1+ei)/ei)*(3*p/hs)^(1-nparam)/(3 + a*a - a*sqrt(3.)*((ei0-ed0)/(ec0-ed0))^alpha);

% fs = 2*(1+e)/(3*r_parms)*(1/exp1 - 1/exp2);
% fs = 2*(1+e)/(3*r_parms)*(1/exp_lambda);
fs = 2*(1+e)/(3*1.0*r_parms)*(1/exp1);
% New definition tensors - Ene 2025
Lep = I_pe*fs*(I*(3*p) + f_v*T*kron_delta + ag^2*(T*T')/(3*p));
Nep = I_pe*fs*fd*(ag*(T + T_star));

% complete von wolferdorf model
% Lep= fs*(F1*F1*I + a*a*T_str*T_str')/stress_trace(TT_str);
% Nep= fs*fd*F1*a*(T_str+T_str_star)/stress_trace(TT_str);

% Wang, Wu, Zhang and Kim (2020)- Overconsilidated clays
% lambda_oed= 0.06;
% kappa     = 0.025;
% yy = (lambda_oed+kappa)/(lambda_oed-kappa);
% c1 = 2*sqrt(3)*ag*(yy-1)/(9*r_parms);
% c2 = 1/3*(sqrt(3)*ag*yy-ag^2-3);
% c3 = ag^2+3*(1-c1);
 
% fs = 3*sqrt(3)*p/(2*ag)*(1/lambda_oed + 1/kappa);
% Lep = (fs)*(I*c1 + c2*T_str*kron_delta + c3*(T_str*T_str'));

%Flow rule m and limit stress condition Y
% qqa = F^2-ag^2*Tstrstarnorm^2/(F^2*(1+c2) + ag^2*sqrt(T_str'*M*T_str)^2);
% m = -T_str_star - qqa*T_str;

% I1= p*3;
% I2= T(1)*T(2)+ T(2)*T(3) + T(3)*T(1);
% I3= T(1)*T(2)*T(3);
% Y = ((yy-1)*(I1*I2/I3-9)+8*tan(phic)^2)/(8*yy*tan(phic)^2);
% fs = 3*sqrt(3)*p/(2*ag)*(1/lambda_oed - 1/kappa);
% Nep = fs*(*m/sqrt(m'*I*m));

% e   = ec;

% using factors by wu
% Lep = (fs)*(I*(3*p_str) + T_str*kron_delta + ag^2*T_str*T_str'/(3*p_str));
% Nep = (fs)*fd*(ag*(T_str + T_str_star));
% e   = ec0*exp(expe);

% using Von Wolferdorf f_s
% Lep = (fs)*I_pe*(C1*I*(3*p_str) + f_v*C2*T_str*kron_delta + C3*T_str*T_str'/(3*p_str));
% Nep = (fs)*I_pe*fd*(C4*(T_str + T_str_star));

% complete wu bahuer
% Lep= fs/I_pe*(a*a*I + T_str*T_str')/stress_trace(TT_str);
% Nep= fs/I_pe*fd*a*(T_str+T_str_star)/stress_trace(TT_str);

global log_internal_vars log_index
log_internal_vars.fs(log_index) = fs;
log_internal_vars.I_pe(log_index) = I_pe;
log_internal_vars.f_v(log_index) = f_v;
log_internal_vars.a01(log_index) = a01;
log_internal_vars.a0(log_index) = a0;
log_internal_vars.ag(log_index) = ag;

%% Post liquefaction parameters
cz          = parms(21);
zmax        = parms(22);
cl          = parms(23);
nl          = parms(24);
pth         = parms(25);
cr          = parms(26);
lambda1     = parms(27);
lambda2     = parms(28);
% lambda3     = parms(32);

Mc = 6*sin(phic)/(3-sin(phic));
Fd = (qdev/p)/(Mc*fd/a01)-1;
% Fd = (qdev/p)/(Mc*fd)-1;
% Fd = (qdev/p)/(Mc*fd*F1)-1;


% Equations Fabric change effect
if sqrt(deps'*I*deps)  ~=  0
        dev_strain_rate_norm = sqrt(de'*I*de);        
        dev_strain_rate = de/dev_strain_rate_norm;
        norm_deps = sqrt(deps'*I*deps);        
    if Fd > 0 %dilatation zone implies that z envolves only after phase transformation
%         z_dot = 0*Fd*norm_deps*(zmax*dev_strain_rate - zint);
        z_dot = cz*Fd*norm_deps*(zmax*dev_strain_rate - zint);
    else
        z_dot = zeros(6,1);
    end
        
else
    dev_strain_rate_norm = 0;
    dev_strain_rate = zeros(6,1);
    z_dot = zeros(6,1);
end


%% new state variable
if zint'*deps > 0 && sqrt(zint'*zint) >= Zdd
%     dZZ = 0*sqrt(z_dot'*z_dot);
    dZZ = sqrt(z_dot'*z_dot);
%     dZZ = sqrt(zint'*zint);
else
    dZZ = 0;
end

% FF = (1 - dZZ/zmax);
FF = (1 - Zdd/zmax);
%%
if zint'*I*dev_strain_rate >= 0
    Ad = 1;
else
    Ad = 1 + -zint'*I*dev_strain_rate;
% Ad = 1;
end


B  = Lep\Nep;
Bp = 1/3*kron_delta'*stress_trace(B); % I*tr(B)/3
Bs = B-Bp;
Bnew = Ad*Bp + Bs;

if Dv == 0
Nep  =  Lep*Bnew;  %signo negativo en -Lep (Ver artículo de Liao et al., 2022)
end

% Dc = -B;

%Semifluized state
pr = p/(pth);
fl = 0.01;
% fl = 0;

if pr <= 1
    l_dot   = abs(des)*( cl*(1-pr )*(1-lint)^nl ) - cr*lint*abs(dev);
% l_dot   = 0.0*abs(des)*( 0*(1-pr )*(1-lint)^0 ) - 0*lint*abs(dev);
l_dot2  = 0*abs(Nep'*deps)*( 0.3*cl*(1-pr )*(1-lint2)^nl ) - cr*lint*abs(dev);
    
    
lambda= lambda1*exp(lambda2*((ei-e)/(ei-ed)));

%     lambda = lambda1*exp(((ei-e)/(ei-ed)));  %Corregir se fuera necesario
    % bbb = 0.9;
    bbb = 1;

    Se    = (1 - (1-pr))^(lambda*bbb*lint) + fl;
%     Se    = 1;

else
    
    l_dot = 0;
    l_dot2= 0;
    Se    = 1;
end

% if DDs > eq_lim
% Br= En/(En+eq_lim);
if Dv > 0
    Dv = Dv;
else
    Dv = 0;
end

dE= Dv*p/(3*p)+qdev*Ds/(3*p); %particle breakage
d_strain2 = dev_strain_rate_norm;
% d_strain2 = parms(20) *(1 - Strain_total2)*dev_strain_rate_norm;
% e = e*exp(Br*e_br);

% end
% dE = 0;
Hep = zeros(6,1);
Hep(1:3) = -(1+e);




% log_internal_vars.z_dot(log_index) = z_dot;
log_internal_vars.eps_a(log_index) = y(3); 
% log_internal_vars.Se(log_index) = Se;
log_internal_vars.FF(log_index) = FF;
log_internal_vars.dZZ(log_index) = dZZ;
log_internal_vars.qdev(log_index) = qdev;
log_internal_vars.p(log_index) = p;
log_internal_vars.pr(log_index) = pr;
% log_internal_vars.l_dot(log_index) = l_dot;
% log_internal_vars.lint(log_index) = lint;
%also
log_internal_vars.Fd(log_index) = Fd;
log_internal_vars.dev_strain_rate_norm(log_index) = dev_strain_rate_norm;
log_index = log_index + 1;

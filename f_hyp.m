function F = f_hyp(y,S,E,V,parms)

% F_EP: computes RHS of the evolution ODE for plastic processes
%

M=[1,  0,  0,  0,  0,  0
   0,  1,  0,  0,  0,  0
   0,  0,  1,  0,  0,  0
   0,  0,  0,  2,  0,  0
   0,  0,  0,  0,  2,  0
   0,  0,  0,  0,  0,  2];

M2=[1,  0,  0,  0,    0,    0
   0,  1,  0,  0,    0,    0
   0,  0,  1,  0,    0,    0 
   0,  0,  0,  0.5,  0,    0
   0,  0,  0,  0,    0.5,  0
   0,  0,  0,  0,    0,    0.5];
kron_delta = [1;1;1;0;0;0];

mR	   	= parms(12);
ny 		= max(size(y));
qint 		= y(13:ny,1);
istrain		= qint(2:7,1);

ny = max(size(y));
F = zeros(ny,1);

[Lep] = hyp_stif(y,parms,[0,0,0,0,0,0]');
Hdel = zeros(6,6);

Maux=Lep*mR;
A    = S*Maux+E;
depsaux = A\V;

[Lep,Nep,Hep,Ftensor_dot,H, Se, dE,p,FF,dZZ,~,ldot2] = hyp_stif(y,parms,depsaux);

% B = -Lep\Nep;
% norm_deps = sqrt(depsaux'*M2*depsaux);
% norm_B    = sqrt(B'*M*B);
% norm_Bdeps= sqrt((B'*M2*depsaux)^2);
% cos_w= norm_Bdeps/(norm_deps*norm_B);
% % 
% Nep = Nep*cos_w^2;

% % --- Toggle plástico global (ej.: parms(34)=0 apaga)
% plastic_on = (numel(parms) >= 34) && (parms(34) ~= 0);
% if ~plastic_on
%     Nep(:,:) = 0;   % corta la rama “plástica” hipoplástica
% end

%% compute A matrix and solve for deps

Maux=Lep*mR;
A    = S*Maux+E;
depsaux = A\V;

[Mep,Hdel,dStrain,dStrain2]=M_istr(Lep,Nep,depsaux,istrain,parms,y(28,1),FF,y(30,1));
Mep = Mep*(1/3*(1-Se)*(kron_delta.*kron_delta') + Se*diag(ones(6,1)));


A    = S*Mep+E;
depsMNR = A\V;
depsinit = depsMNR;
sigpresc = Mep*depsMNR;

[Mep,Hdel,dStrain,dStrain2]=M_istr(Lep,Nep,depsMNR,istrain,parms,y(28,1),FF,y(30,1));
Mep = Mep*(1/3*(1-Se)*(kron_delta.*kron_delta') + Se*diag(ones(6,1)));

%% Modified Newton-Raphson interations for correction of deps
stresscontrol=zeros(6,1);
for i = 1:6
  for j = 1:6
    if  (S(i,j) > 0)
      stresscontrol(j)=1;
    end
  end
end

isstresscont=sqrt(stresscontrol'*M*stresscontrol);

if isstresscont > 0
  tol_mnriter  = parms(3);
  itererror=1000;
  
  while itererror > tol_mnriter
  
        sigest= Mep*depsMNR;
        unbalsig=sigpresc-sigest;
 
        for i = 1:6
	  if  (stresscontrol(i) == 0)
	    depsMNR(i)=depsinit(i);
	    unbalsig(i)=0;
	  end
        end

        depsMNR=depsMNR + Mep\unbalsig;
        itererror=sqrt(unbalsig'*M*unbalsig);

        [Mep,Hdel,dStrain,dStrain2]=M_istr(Lep,Nep,depsMNR,istrain,parms,y(28,1),FF,y(30,1));
        Mep = Mep*(1/3*(1-Se)*(kron_delta.*kron_delta') + Se*diag(ones(6,1)));

  end
end
%% F vector

F(1:6,1)   = depsMNR;
F(7:12,1)  = Mep*depsMNR;
de= Hep'*depsMNR;
distrain=Hdel'*depsMNR;
F(13,1) = de;
F(14:19,1) = distrain;
F(20:25,1) = Ftensor_dot;
F(26,1) = H;
F(27,1) = dE;
F(28,1) = dStrain;
F(29,1) = dZZ;
F(30,1) = dStrain2;
F(31,1) = ldot2;




function [SS,EE,INV_S,INV_E,HARD] = updateModel(y0,parms,nspb,path_info)

% UPDATE: compute the evolution of strain, stress and internal variables
%         over the prescribed loading path
%
% version 1.0 for rate-independent elastoplastic models
%
% integration methods implemented:
% explicit adaptive Runge-Kutta (3rd order) with error control 
%
% written by C. Tamagnini - jan. 2012

%% initialize load rate vector

V=zeros(6,1);

%% recover input data

tol_f=parms(1);

ny = max(size(y0));

eps  = y0(1:6,1);
sig  = y0(7:12,1);
qint = y0(13:ny,1);

icode = path_info(1,:);
nstep = path_info(2,:);
DX    = path_info(3,:);


dx = DX./nstep;  % increment size for each sp branch

%% copy initial state into global output arrays

kstep = 1;

[p,q,z] = inv_s(sig);
[epsv,epss] = inv_e(eps);

SS(kstep,:)=sig';
EE(kstep,:)=eps';
INV_S(kstep,:)=[p,q,z];
INV_E(kstep,:)=[epsv,epss];
HARD(kstep,:)=[qint'];

%% loop over the SPB

for i=1:nspb
    
    disp(' ')
	disp(['STRESS PATH BRANCH # ', int2str(i)])
	disp(' ')

	n=nstep(i);
	k=icode(i);
    dload = dx(i);

% loop over load steps

	for j=1:n

		kstep=kstep+1;
		disp(['STEP # ',int2str(kstep-1)])
 
% load vector V

        V(6)=dload;

% initial state

    y_k=[eps',sig',qint']';

%  hypoplastic evolution equation

        y_k = hypo_update(y_k,V,parms,k);

% recover state variables and invariants at the end of substepping

		eps  = y_k(1:6,1);
		sig  = y_k(7:12,1);
		qint = y_k(13:ny,1);

        [p,q,z] = inv_s(sig);
        [epsv,epss] = inv_e(eps);
        
% move final states into stress and strain path vectors
%         plot(eps(1),sig(1)-sig(3),'o'); hold on; pause(0.0002)
        SS(kstep,:)=sig';
        EE(kstep,:)=eps';
        INV_S(kstep,:)=[p,q,z];
        INV_E(kstep,:)=[epsv,epss];
        HARD(kstep,:)=[qint'];

% bottom of the loop over load steps

	end

% bottom of the loop over stress path branches
 
end

% INPUT_DATA: input material constants and store them in vector parms
%             set loading path info
%
% nspb          = no. of stress path branches;
% icode(1:nspb) = vector of path type codes (one for each branch)
% nstep(1:nspb) = vector of step numbers (one for each branch)
% DX(1:nspb)    = vector of total "load" increment for each branch
%
% loading path codes
%
% 1. strain controlled undrained TX compression (axis direction: x_3)
% 2. strain controlled ED compression (axis direction: x_3)
% 3. stress controlled drained TX compression (axis direction: x_3)
% 4. mixed control drained TX compression
% 5. mixed control ED compression
% 6. mixed control plane strain compression
% 7. strain controlled undrained simple shear
% 8. stress (mixed) controlled simple shear

%% material parameters
% 1.0486 0.1075 36.0093 0.8325 1.1142 0.5526
% phic     	= 32.6328;
phic     	= 37;
hs       	= 3.4156e5;
hs_reverse	= 1.0;
nparam     	= 0.4564;
ed0 		= 0.55; % 0.54
ec0	     	= 1.0345;
ei0	   	    = 0.9775;
alpha	   	= 0.62;
beta	   	= 4.875;
r_parm      = 1.08;
c           = 41*0;



% Liquefaction params

cz      = 1500;
zmax    = 15;
cl      = 50;
nl      = 8;
pth     = 10*1;
cr      = 1000;
lambda1 = 0.15;
lambda2 = 3.0;


q_br = 0.5;
kb_br= 0.0128;
a_br = 15;

mR = 2.5;
mT = 2.5;
Rparam = 5.0e-5*2; % *8
betar = 0.15; % 0.18
chi = 5.5*1.10; % 3

C3 = 1.1647;


max_ksub 	= 1000;       % max no. of substeps allowed
err_tol  	= 1.0e-4;     % normalized error tolerance for adaptive integration  
tol_mnriter	= 1.0e-10;     % tolerance on Modified Newton Raphson estimation of deps

parms=[max_ksub;
       err_tol;
       tol_mnriter;
       phic;
       hs;
       nparam;
       ed0;
       ec0;
       ei0;
       alpha;
       beta;
       mR;
       mT;
       Rparam;
       betar;
       chi;
       hs_reverse;
       c;
       r_parm;
       q_br;
       [cz zmax cl nl pth cr lambda1 lambda2]';
       kb_br;
       a_br;
       C3];

   
   
   
%% type and characteristics of loading path

nspb=1;

% icode=[5,5,5];
icode=[4]*ones(1,nspb);%[3,3,3,3,3];
nstep=[500]*ones(1,nspb)%[500,500,500,500,500];
% DX=[787.5,-750,750];
% DX=[30,-30*2,30*2];
DX= 0.25*[1 2*ones(1,nspb-1).*(-1).^(1:nspb-1)];
% DX= [800 -762.5 762.5];

path_info = [icode;
             nstep
             DX];

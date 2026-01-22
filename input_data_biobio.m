
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
%% ORIGINAL HYPOPLASTIC MODEL

%Le pegan bien tcu
phic     	= 34.1;
hs       	= 21.18e6;
hs_reverse	= 1.0;
nparam     	= 0.23
ed0 		= 0.610;
ec0	     	= 0.947;
ei0	   	    = 1.09;

alpha	   	= 0.95%0.97%0.7;
r_parm      = 1.13%1.7;
c           = 0.0; 
beta	   	= 0.13%0.01

mR = 1.9;
mT = 1.41%0.77*mR;
Rparam = 1.0e-4*1.0;
betar = 0.1; % 
chi = 3.0;%22; 
%% LIQUEFACTION PARAMETERS
init_state_biobio
%Fabric change 
cz      = 700;
zmax    = 20;
%Default constants model
nl      = 8;
pth     = 10.0;
%Semifluized state
cl      = 60;
cr      = 1000;
lambda1 = 0.2;
lambda2 = 3.0;

%% NEW STATE VARIABLES
chimax= 20;
eta=0.1

q_br = 2.0%14.5;%9.9575; % parms(20)= beta2 
kb_br=  0.0740;
a_br = 9.635; %parms(30) = eps_acc1

%Y estos?
C3 = 2.0;
C5 = 0.08;
%% PARAMS DEFINITION
max_ksub 	= 1000;       % max no. of substeps allowed
err_tol  	= 1.0e-4;     % normalized error tolerance for adaptive integration  
tol_mnriter	= 1.0e-10;    % tolerance on Modified Newton Raphson estimation of deps

parms=[max_ksub;%1
       err_tol;%2
       tol_mnriter;%3
       phic;%4
       hs;%5
       nparam;%6
       ed0;%7
       ec0;%8
       ei0;%9
       alpha;%10
       beta;%11
       mR;%12
       mT;%13
       Rparam;%14
       betar;%15
       chi;%16
       hs_reverse;%17
       c;%18
       r_parm;%19
       q_br;%20
       %21 %22 %23 %24 %25 %26  %27    %28
       [cz zmax cl  nl pth  cr lambda1 lambda2]';
       kb_br;%29
       a_br;%30
       C5;%31
       e0;%32
       chimax;%33
       eta;
       ]
%% TYPE AND CHARACTERISTICS OF LOADING
nspb = 30%58%58
% icode=[5,5,5];
icode=[1]*ones(1,nspb);%[3,3,3,3,3];
nstep=[500]*ones(1,nspb);%[500,500,500,500,500];

% DX=[700 -500 700];
DX= 0.0005*[1 2*ones(1,nspb-1).*(-1).^(1:nspb-1)];
% DX= 0.3*[1 2*ones(1,nspb-1).*(-1).^(1:nspb-1)];
path_info = [icode;
             nstep;
             DX];

 
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

% phic     	= 34;
% hs       	= 3.4156e5%  3.97407e8 %3.4156e8 %[100-500]MPa
% hs_reverse	= 1.0;
% nparam     	= 0.233 %[0.2 - 0.4]. 
% ei0	   	    = 0.941 %1.09 %; % Herle and Gudehus (1999) --> eio = 1.15ecr; eio=1.2eco
% ec0	     	= 0.824 %0.948 %ei0/1.2   %ei0/1.2; %0.948  
% ed0 		= 0.508 %0.610  
% alpha	   	= -1.94 %0.831 %0.274  %  [0.05 - 0.15].
% beta	   	= 0.53 %0.5 %1.46 [0.5 - 1.5].
% r_parm      = 1.5%2.15 
% c           = 0
% % intergranular strain
% mR = 4.6*.5;
% mT = 0.55*mR;
% Rparam = 5.e-5*2; % *8
% betar = 0.15; % 0.18
% chi = 1.5; % 3
%% Final parameters karlsuhe sand
 %     hs        n    alpha    beta  r_parm    Rparam       q_br     
% x0 = [4.0e6   0.27   0.8311     0.11    1.01    1.0e-04   9.95]
phic     	= 32.5%32.5%32.5;
hs       	= 4e6 
hs_reverse  = 0.05;
nparam     	= 0.27;
ed0 		= 0.677; % 0.54
ec0	     	= 1.054;
ei0	   	    = 1.212%*0.875;
% ei0	   	    = 1.0908;
alpha	   	= 0.8311
beta	   	= 0.1114 %0.85 para original model response 
r_parm      = 1.0902;
c           = 41*0;

mR = 2.2;
mT = 1.1;
Rparam = 1.0e-4; 
betar = 0.15;
chi = 5.0; 
%% Liquefaction params
init_state_karlsuhe
Dr      = (ei0-e0)/(ei0-ed0);
% Xinew   = -5.421 + 8.231*Dr + 0.0447*sum(sig0)/3;
% Xinew   = 10
% lambda1 = 1.4;
% lambda2 = 4.349*(sum(sig0)/3)^-0.3924 + 0.3361;% anisotropic a01 =1
% lambda3 = -1.231 + 2.107*Dr + 0.03207*(sum(sig0)/3) - 0.03016*Dr*(sum(sig0)/3) -2.644e-05*(sum(sig0)/3)^2;
%% Fabric change effect
cz = 500;
% zmax    = min(20*Dr^3.85,15);%  1 < Dr 66 ;(10 Dr 94)
% zmax    = min(30.4*Dr^6.31,15);%

zmax = 99.68*(Dr^-2.606)*(0.25^2.451);
% zmax = 25.3;
%% Semifluized state
%Present study
cl      = 50; % [10:66 7:77 5:94 3:98]
cr      = 1000;
% lambda1 = 7.495*exp(-2.474*Dr)*(0.25^1.614);
lambda1 = 0.2;
lambda2 = 3.0;
%% NEW STATE VARIABLES
% chimax = 18.1
chimax = 59*(Dr^-2.064)*(0.25^1.506);
q_br =  9.9575; %* parms(20) = beta2 (monotonic regimen)
%%
kb_br=  0.0740;
a_br = 9.635; 
C3 = 2.0;
C5 = 0.08;
%% Total parameters
%Default parameters
nl      = 8;
pth     = 10;
%%
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
       %21  22  23 24 25  26  27      28
       [cz zmax cl nl pth cr lambda1 lambda2 ]';
       kb_br;%29
       a_br;%30
       C5;%31
       e0;%32
       %Xinew;
       chimax
       ];      
%% type and characteristics of loading path

nspb = 20%181%55
% icode=[5,5,5];
icode=[1]*ones(1,nspb);%[3,3,3,3,3];
nstep=[500]*ones(1,nspb)%[500,500,500,500,500];

DX= 0.0005*[1 2*ones(1,nspb-1).*(-1).^(1:nspb-1)]; %TCU
% DX= 0.3*[1 2*ones(1,nspb-1).*(-1).^(1:nspb-1)]; %TCU


% DX= [700 -700 ];

path_info = [icode;
             nstep
             DX];
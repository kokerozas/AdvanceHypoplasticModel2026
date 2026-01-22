%% Driver to test liquefaction model and modification of I_d density factor
clear all 
clc

input_data_karlsuhe;
init_state_karlsuhe;

% input_data_biobio;
% init_state_biobio;

init_data = load('D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2023\TMU - triaxial-monotonic-undrained\TMU_fine_sand.txt');

kk = 1;

for i = 3:5:length(init_data)

path_info(3,1) = init_data(i,1);  %load increment
y0(7:9,1)      = init_data(i,2);  %sigmas
y0(13,1)       = init_data(i,3);  %void ratios

% parms(10)      = 0.6;
% parms(11)      = 1.5; b

[SS,EE,INV_S,INV_E,HARD] = updateModel(y0,parms,nspb,path_info);

figure(8)
xyplot_vmh; hold on;

name = ['D:\Driver_hypo_sand\hypo_sand_istrain_new_lique2023\TMU - triaxial-monotonic-undrained/TMU' num2str(i) '.dat'];

a = importdata(name);

eps1 = a.data(:,1)/100;
p = a.data(:,7);
q = a.data(:,6);
epsv = a.data(:,2)/100;
e = a.data(:,5);

subplot(2,2,1)
plot(eps1,q,'b-'); hold on;
subplot(2,2,2)
plot(p,q,'b-'); hold on;
subplot(2,2,3)
plot(eps1,epsv,'b-'); hold on;
subplot(2,2,4)
plot(epsv,q,'b-')

[~, pos] = unique(eps1);

epsv_in(:,kk) = interp1(eps1(pos),epsv(pos), EE(:,3));
q_in(:,kk)    = interp1(eps1(pos),q(pos), EE(:,3));

kk =  1+ kk;
end
% 
% writematrix(epsv_in,'eps_in_karlsuhe.txt')
% writematrix(q_in,'q_in_karlsuhe.txt')
%%
TAU    = load('q_in_karlsuhe.txt');
DEV    = load('eps_in_karlsuhe.txt');
% plot(eps1,q,'b');hold on

plot(DEV,TAU,'b');hold on
% plot(INV_E(:,1),INV_S(:,2),'r')
%% 

input_data_karlsuhe;
init_state_karlsuhe;

% TAU    = load('q_300.txt');
% DEV    = load('epsv_300.txt')

ObjFun = @(x) fitnessfunction_karlsuhe(x,TAU,DEV,init_data);

% options = optimset('TolFun',0.8,'TolX',1,'OutputFcn', @myoutput);
options = optimset('OutputFcn', @myoutput);
options.MaxFunEvals = 1000;
options.MaxIter     = 1000;
options.UseParallel = 1;

%                             beta
x0 = [0.9234*0.9    1.0902    0.1114    9.9575];
% x0 = [0.696    7.9604*1    1.0491*0.7    0.6226/2    0.0740    9.6350    2 0.08]
% 0.7630    6.3015    1.0137    0.3303    0.0778    8.8876    1.9587    0.0877    1.0749

% alpha(ee)  = -0.482*x^18.63+0.6103;
% beta(ee)   = 1.25*x^16.63+6.344;
% rparam(ee) = (1.531*x^2 + -2.366*x + 0.9191) / (x^2 + -1.565*x + 0.6181);
% q_br(ee)   = 0.609*exp(-0.8411*x) + 1.278e-16*exp(35.28*x)
% En(ee)     =  0.1864 + 0.01375*cos(x*19.63) + -0.1*sin(x*19.63)

% x0 = [1.1757    5.1615    1.2428   19.3490    0.0979   16.8263    0.8719]
x1 = [0.5 0  0.5 0.1  0     10  ];
x2 = [1.2 12  1.7 3.0  0.15  35  ];
% A  = [0 0 0 0 0 0  1 -1  0;
%       0 0 0 0 0 0 -1  0  1];
% b  = [0;0];
[x , fval] = fminsearch(ObjFun,x0,options);
% [x , fval] = patternsearch(ObjFun,x0,[],[],[],[],x1,x2,[],options);
%% 

input_data_karlsuhe;
init_state_karlsuhe;

[SS,EE,INV_S,INV_E,HARD] = updateModel(y0,parms,nspb,path_info);
% [SS,EE,INV_S,INV_E,HARD] = updateModel_stress(y0,parms,nspb,path_info,35,35);

% semilogx(SS(:,3),HARD(:,1)); hold on;

figure(3)
subplot(1,2,1); hold on;
plot(EE(:,3),SS(:,3)-SS(:,1),'b');
xlabel('\epsilon_{s}'); ylabel('q (kPa)')
% xticks([-0.005:0.0025:0.005]*2)
% plot(EE(:,6),SS(:,6),'b');
subplot(1,2,2); hold on;
plot(INV_S(:,1),SS(:,3)-SS(:,1),'b');
xlabel('p (kPa)'); ylabel('q (kPa)')
xticks([0:25:100])




xyplot_vmh; hold on;
% 
%%


% 
% Linear model Poly2:
%      f(x) = p1*x^2 + p2*x + p3
% Coefficients (with 95% confidence bounds):
%        p1 =      -4.487  (-10.64, 1.67)
%        p2 =       5.728  (-0.6708, 12.13)
%        p3 =     -0.2381  (-1.842, 1.365)
% 
% Linear model Poly2:
%      f(x) = p1*x^2 + p2*x + p3
% Coefficients:
%        p1 =      -5.052
%        p2 =       4.326
%        p3 =       1.173
%        
% Linear model Poly2:
%      f(x) = p1*x^2 + p2*x + p3
% Coefficients (with 95% confidence bounds):
%        p1 =      -2.587  (-5.216, 0.04171)
%        p2 =       3.315  (0.6276, 6.002)
%        p3 =     -0.3228  (-0.895, 0.2495)      
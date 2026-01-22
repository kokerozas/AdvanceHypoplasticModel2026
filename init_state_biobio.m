    % INIT_STATE: initialize state variables
%
% eps0  = initial strain state
% sig0  = initial stress state
% qint0 = initial value of internal variables (stored in a 1d array)

%% preallocate arrays

eps0 = zeros(6,1);
sig0 = zeros(6,1);

%% input initial strain

% eps0 = [0.0; 0.0; 0.001; 0.0; 0.0; 0.0]   % uncomment this line if eps0 =/= 0

%% input initial stress
% erel = [0.5876 0.6494 0.7207 0.9347];
% erel = [0.539 0.609 0.675 0.770];
fg   = 2.1;
k0   = 1-sind(phic);
% sig0 = [fg*k0; fg*k0; fg; 0.0; 0.0; 0.0];
sig0 = 1*[100; 100; 100; 0.0; 0.0; 0.0];
e0 	   = ei0 - .77*(ei0-ed0) ;
% e0 	   = ec0 - .25*(ec0-ed0) ;

% e0 = 0.866 %(Para un dr = 21% --> dr=24 con fórmula)  
% e0 = 0.759
istrain0 = -0.8*Rparam/sqrt(3)*[1.0; 1.0; 1.0; 0.0; 0.0; 0.0];
%% input initial internal variables
A = 0;
F0= 0.5;
ldot = 0;
delta= 45;
Mmatrix = [sind(delta)      0        cosd(delta) ; 
            0               1           0; 
            cosd(delta)     0       -sind(delta)];
        
% Fmatrix = diag([ -1/sqrt(6) -1/sqrt(6) 2/sqrt(6)]);
% Fmatrix = diag([ 2/sqrt(6) -1/sqrt(6) -1/sqrt(6)]);
% Ftrans  = Fmatrix*Mmatrix;

% Ftensor = -F0*[ 2/sqrt(6) -1/sqrt(6) -1/sqrt(6)  0 0 0]';
Ftensor = F0*[ -1/sqrt(6) -1/sqrt(6) 2/sqrt(6)  0 0 0]';

% Ftensor = [0 0 0 0 0 0]';
        %                          
qint0 = [e0;istrain0;Ftensor;ldot;0;0;0;0;0];
%% collect all state variables in vector y0

y0 = [eps0',sig0',qint0']';

Relative_density   = 0.90;
ymin               = 1.450;
ymax 			   = 1.771;
yrelative          = ymin*ymax/(ymax-Relative_density*(ymax-ymin));
e_max              = 1.09 %0.947 %0.91;
e_min              = 0.610  %0.51;
e_relative 		   = e_max - Relative_density*(e_max-e_min) ;

% e_rel = 0.7464 60%
% e_rel = 0.5658 95%
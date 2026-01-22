function [p,q,z] = inv_s(sigma)

% INV_S: Invariants of the stress tensor sigma
%        stored as a (6x1) vector in Voigt notation 
%
%        p=tr(sig)/3
%        q=sqrt((3/2)*tr(s*s))
%        z = sin(3*theta), theta = Lode angle
%

fact=1/3;
m=[1;  1;  1;  0;  0;  0];
M=[1,  0,  0,  0,  0,  0
   0,  1,  0,  0,  0,  0
   0,  0,  1,  0,  0,  0
   0,  0,  0,  2,  0,  0
   0,  0,  0,  0,  2,  0
   0,  0,  0,  0,  0,  2];

p=fact*m'*sigma;
if p < 0.0001
    p = 0.0001;
end
s=sigma-p*m;
sn=s'*(M*s);

q=sqrt(3*sn/2);

% components of s*s as a (6x1) vector t1

t1=[s(1)*s(1)+s(4)*s(4)+s(6)*s(6);
    s(4)*s(4)+s(2)*s(2)+s(5)*s(5);
    s(6)*s(6)+s(5)*s(5)+s(3)*s(3);
    2*(s(1)*s(4)+s(4)*s(2)+s(6)*s(5));
    2*(s(4)*s(6)+s(2)*s(5)+s(5)*s(3));
    2*(s(6)*s(1)+s(5)*s(4)+s(3)*s(6))];

% Lode angle

if q==0.0
    z = 0;
else
    temp1=fact*t1'*s;
    z=27.0*temp1/(2.0*q*q*q);
    if abs(z) > 1.0
        z=z/abs(z);
    end
end


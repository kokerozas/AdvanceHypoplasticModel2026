clear all
clc


hs = 4.62e7;
nparam = [0.245 0.35 0.4];
r_params = 0.78;


e = [linspace(0.610,1.09,100)'];
p = [linspace(0,100,100)']


for k = 1:3
    for i = 1:length(e)
        for j=1:length(p)
            exp1 = nparam(k)*e(i)*(3*p(j)/hs).^nparam(k);
            fs{k}(i,j) = 2*(1+e(i))/(3*r_params)*(1/exp1);
        end
    end
end

 [X,Y]= meshgrid(e,p)
 for i=1:1
    surf(X,Y,fs{i}) ; 
    set(gca,'YScale','log')
 end

% plot3(e,r_params,fs,'ro') ; set(gca, 'xscale','log')
    
 xlabel('void ratio, e(-)')
 ylabel('mean pressure, p (kPa)')
 zlabel('fs factor (-)')
   
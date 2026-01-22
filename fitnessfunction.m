function y = fitnessfunction(x,TAU,DEV)
% input_data_karlsuhe;
% init_state_karlshue;

input_data;
init_state;

rel  = [0.17 0.62 0.77 0.90];

    for ik = 1:4
%                 ed0 		= x(9); % 0.54
                ed0 		= 0.55;
                ei0	   	    = 0.9775;
%                 ei0	   	    = 1.05;
%                 ei0	   	    = x(7);
%                 ei0	   	    = 1.02;

%               y0(13)     = ei0 - rel(ik)*(ei0-ed0);
                 y0(13)     = x(4) - rel(ik)*(x(4)-x(3));
                parms(10)  = x(1); %alpha
                parms(11)  = x(2); %beta
                parms(19)  = x(3); %rparms
                parms(20)  = x(4);
                parms(29)  = x(5);
                parms(30)  = x(6);
                
                              
     %comentar línea 29  y 30          
%                 parms(8)   = x(7);
%                 parms(31) = x(8);

          [~,Es,INV_S,INV_E,~] = updateModel(y0,parms,nspb,path_info);
                SS(:,ik) = INV_S(:,2);
                EE(:,ik) = INV_E(:,1);
                
    end

    SS(isnan(SS)) = 0;
    EE(isnan(EE)) = 0;
    
    set_a = 1:400;
    set_b = 1:length(TAU);
    set_b = setdiff(set_b,set_a);
%     pause()
    y1 =  norm(TAU(set_a,:) - SS(set_a,:))/norm(TAU(set_a,:));
    y_v=  norm(DEV(set_a,:) - EE(set_a,:))/norm(DEV(set_a,:));
    y2 =  norm(TAU(set_b,:) - SS(set_b,:))/norm(TAU(set_b,:));
    y_v2= norm(DEV(set_b,:) - EE(set_b,:))/norm(DEV(set_b,:));
    y  = 0.7*(0.5*y1 + 0.5*y_v) + 0.3*(0.5*y2 + 0.5*y_v2)
    
    subplot(1,2,1)
    plot(Es(:,3) , TAU, 'b');hold on
    plot(Es(:,3) , SS , 'r'); hold off
    subplot(1,2,2)
    plot(Es(:,3) , DEV, 'b');hold on
    plot(Es(:,3) , EE,  'r'); hold off
    
    pause(0.1)
%     clf 
end

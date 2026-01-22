function y = fitnessfunction_biobio(x,TAU,DEV,init_data)
input_data_biobio;
init_state_biobio;

indexx = 1:4%5:10:length(init_data)
    for ik = indexx
%                 path_info(3,1) = init_data(ik,1);
%                 y0(7:9,1)      = init_data(ik,2);
                y0(13,1)       = init_data(ik,3);
%               y0(13)         = rel(ik);

%    alpha   r_parm  beta    epsmon
% x0 = [0.77   1.7     0.13     15.3]

                   parms(10)  = x(1); %alpha
                   parms(11)  = x(2); %beta
                   parms(19)  = x(3); %r_param
                   parms(20)  = x(4); %epsmon        

                   
                [Se,Es,INV_S,INV_E,hard] = updateModel(y0,parms,nspb,path_info);
                SS(:,ik) = Se(:,3) - Se(:,1);  
                EE(:,ik) = INV_E(:,1);
                EE3(:,ik) = Es(:,3);
                         
    end

 SS(isnan(SS)) = 0;
 EE(isnan(EE)) = 0;
    
    set_a = 1:470
    set_b = 1:length(TAU)-1;
    set_b = setdiff(set_b,set_a);
 
% ERROR FOR TMD TEST
    y1 =  norm(TAU(set_a,indexx) - SS(set_a,indexx))/norm(TAU(set_a,indexx));
    y_v=  norm(DEV(set_a,indexx) - EE(set_a,indexx))/norm(DEV(set_a,indexx));
    
    y2 =  norm(TAU(set_b,indexx) - SS(set_b,indexx))/norm(TAU(set_b,indexx));
    y_v2= norm(DEV(set_b,indexx) - EE(set_b,indexx))/norm(DEV(set_b,indexx));
    
    y  = 0.7*(0.6*y1 + 0.4*y_v)+ 0.3*(0.6*y2 + 0.2*y_v2)
%     y  = 0.8*(0.9*y1 + 0.1*y_v)+ 0.2*(0.9*y2 + 0.1*y_v2)
%     y  = 0.9*(0.5*y1 + 0.5*y_v)+ 0.1*(0.5*y2 + 0.5*y_v2)

%     subplot(1,2,1)
%     plot(EE3(set_a,indexx) , TAU(set_a,indexx), '-b','linewidth',1.3);hold on
%     plot(EE3(set_a,indexx) , SS(set_a,indexx) , '-r','linewidth',1.3); hold off
%     subplot(1,2,2)
%     plot(EE3(set_a,indexx) , DEV(set_a,indexx), '-b','linewidth',1.3);hold on
%     plot(EE3(set_a,indexx) , EE(set_a,indexx),  '-r','linewidth',1.3); hold off
 
    set_a = 1:500;
    figure(3)
    subplot(1,2,1)
    plot(EE3(set_a,indexx) , TAU(set_a,indexx), 'b');hold on
    plot(EE3(set_a,indexx) , SS(set_a,indexx) , 'r'); hold off
    subplot(1,2,2)
    plot(EE3(set_a,indexx) , DEV(set_a,indexx), 'b');hold on
    plot(EE3(set_a,indexx) , EE(set_a,indexx),  'r'); hold off
%     
    pause(0.1)
%     clf
      
end

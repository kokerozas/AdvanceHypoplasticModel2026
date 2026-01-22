function y = fitnessfunction_biobio(x,TAU,DEV,init_data)
input_data_biobio;
init_state_biobio;

indexx = 1:4
    for ik = indexx
                path_info(3,1) = init_data(ik,1);
                y0(7:9,1)      = init_data(ik,2);
                y0(13,1)       = init_data(ik,3);
%               y0(13)         = rel(ik);

      % phic % hs     n    alpha  beta  r_param   c     mR    mT     Rparam  betar  chi
%x0 = [ 34 3.96e7  0.233   1.4   0.3    1.5    5.0  2.53   0.77*mR  1.e-4  0.25 1.35];

     %phi    hs     n      alpha  beta   r_param   c      mR     mT   Rparam   betar   chi    q_br
% x0 = [34   4.628e7  0.25   0.7     0.2   1.2     2.6    3.5    2.5    1e-4     0.7    2.5    10];

                parms(4) = x(1)%phi
                parms(5)  = x(2); %hs
                parms(6)  = x(3); %nparm
                parms(10)  = x(4); %alpha
                parms(11)  = x(5); %beta
                parms(19)  = x(6); %r_param
                parms(18)  = x(7); %c
                parms(12)  = x(8); %mR
                parms(13)  = x(9); %mT
                parms(14)  = x(10); %Rparam       
                parms(15)  = x(11); %betar
                parms(16)  = x(12); %chi
                parms(20)  = x(13); %q_br
  
                
                [Se,Es,INV_S,INV_E,hard] = updateModel(y0,parms,nspb,path_info);
                SS(:,ik) = Se(:,3) - Se(:,1);  
                EE(:,ik) = INV_E(:,1);
                EE3(:,ik) = Es(:,3);
                         
    end

 SS(isnan(SS)) = 0;
 EE(isnan(EE)) = 0;
    
    set_a = 1:490
    set_b = 1:length(TAU)-1;
    set_b = setdiff(set_b,set_a);

 
% ERROR FOR TMD TEST
    y1 =  norm(TAU(set_a,indexx) - SS(set_a,indexx))/norm(TAU(set_a,indexx));
    y_v=  norm(DEV(set_a,indexx) - EE(set_a,indexx))/norm(DEV(set_a,indexx));
    y2 =  norm(TAU(set_b,indexx) - SS(set_b,indexx))/norm(TAU(set_b,indexx));
    y_v2= norm(DEV(set_b,indexx) - EE(set_b,indexx))/norm(DEV(set_b,indexx));
    y  = 0.8*(0.6*y1 + 0.4*y_v)+ 0.2*(0.6*y2 + 0.4*y_v2)
%         y  = 0.7*(0.5*y1 + 0.5*y_v)+ 0.3*(0.5*y2 + 0.5*y_v2)
    figure(1)
    subplot(1,2,1)
    plot(EE3(set_a,indexx) , TAU(set_a,indexx), '-b','linewidth',1.3);hold on
    plot(EE3(set_a,indexx) , SS(set_a,indexx) , '-r','linewidth',1.3); hold off
    subplot(1,2,2)
    plot(EE3(set_a,indexx) , DEV(set_a,indexx), '-b','linewidth',1.3);hold on
    plot(EE3(set_a,indexx) , EE(set_a,indexx),  '-r','linewidth',1.3); hold off
 

%     figure(3)
%     subplot(1,2,1)
%     plot(EE3(set_a,indexx) , TAU(set_a,indexx), 'b');hold on
%     plot(EE3(set_a,indexx) , SS(set_a,indexx) , 'r'); hold off
%     subplot(1,2,2)
%     plot(EE3(set_a,indexx) , DEV(set_a,indexx), 'b');hold on
%     plot(EE3(set_a,indexx) , EE(set_a,indexx),  'r'); hold off
%     
%     pause(0.1)
%     clf
      
end

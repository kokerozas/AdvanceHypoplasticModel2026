function y = fitnessfunction_karlsuhe(x,TAU,DEV,init_data)
input_data_karlsuhe;
init_state_karlsuhe;

% Lvar = 4:5:length(init_data);

% rel = [0.97 0.8476 0.8136 0.7341 0.697];
indexx = 1:3%3:5:length(init_data);
% indexx = indexx(5);
% indexx = [2 4 7 9 12 14 17 19 22 24];

    for ik = indexx

                path_info(3,1) = init_data(ik,1);
                y0(7:9,1)      = init_data(ik,2);
                y0(13,1)       = init_data(ik,3);
%                 y0(13)         = rel(ik);
                

               
                parms(5)   = x(1); %hs
                parms(6)   = x(2); %nparm
                parms(10)  = x(3); %alpha
                parms(11)  = x(4); %beta
                parms(19)  = x(5); %r_param

                parms(12)  = x(7); %mR
                parms(13)  = x(8); %mT
                parms(14)  = x(9); %Rparam
                parms(15)  = x(10); %betar
                parms(16)  = x(11); %chi
                parms(20)  = x(4); %9.9575

                
              [Se,Es,INV_S,INV_E,~] = updateModel(y0,parms,nspb,path_info);
              SS(:,ik) = Se(:,3)-Se(:,1);
              EE(:,ik) = INV_E(:,1);
              EE3(:,ik) = Es(:,3);
                
    end

    SS(isnan(SS)) = 0;
    EE(isnan(EE)) = 0;
    
    set_a = 1:430;
    set_b = 1:length(TAU)-1;
    set_b = setdiff(set_b,set_a);
%     pause()


    y1 =  norm(TAU(set_a,indexx) - SS(set_a,indexx))/norm(TAU(set_a,indexx));
    y_v=  norm(DEV(set_a,indexx) - EE(set_a,indexx))/norm(DEV(set_a,indexx));
    y2 =  norm(TAU(set_b,indexx) - SS(set_b,indexx))/norm(TAU(set_b,indexx));
    y_v2= norm(DEV(set_b,indexx) - EE(set_b,indexx))/norm(DEV(set_b,indexx));
    y  = 0.7*(0.6*y1 + 0.4*y_v)+ 0.3*(0.6*y2 + 0.4*y_v2)
    
    
    set_a = 1:500;
    figure(3)
    subplot(1,2,1)
    plot(EE3(set_a,indexx) , TAU(set_a,indexx), 'b');hold on
    plot(EE3(set_a,indexx) , SS(set_a,indexx) , 'r'); hold off
    subplot(1,2,2)
    plot(EE3(set_a,indexx) , DEV(set_a,indexx), 'b');hold on
    plot(EE3(set_a,indexx) , EE(set_a,indexx),  'r'); hold off
    
    pause(0.1)
%     clf


      
end

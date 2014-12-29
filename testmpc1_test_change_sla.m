classdef testmpc1_test_change_sla
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
              
        function test_change_sla(o)
            % o.setup_hosts();
            o.num_host=4;
            o.cap=[8 8 8 8]';
            o.speed_ratio=[1 1 1 1]';
            RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',sum(10)));
            o.cost=round(rand(o.num_host,1)*10)+1;
            
            % before change, steady state , initial state
            % after change
            o.num_classes = 1;
            o.num_service=1;
            call=zeros(o.num_service+o.num_classes,o.num_service+o.num_classes);
            call=[0 1
                0 0];
            o.c=rand(o.num_host,o.num_service)+10;
            o.c= [10.49851
                10.2248
                10.1981
                10.7605];
            o.call_graph = call;
            Y=calculate_contact_level(o,call);
            D =  [0.04 ];
            o.d= Y .* (ones(o.num_classes,1) * D);
            T=15;
            
            % before the service arrives
            N=[250];
            Z=[1 ]';
            RT_sla=[0.4 ]';
            x0=calculate_throuput_objective(o, N(:,1), Z, RT_sla, 1);
            [theta,thetaV]=o.solve_nfm(x0);
            thetaV_0 =thetaV; % zeros(o.num_host , o.num_service);
            thetaV
            %  print_lat_mat(o, thetaV)
            
            % over time
            N=[ones(1,T)*250];
            RT_sla=[0.14 ]';
            X_sla = calculate_throuput_objective(o, N, Z, RT_sla, T);
            X_target = X_sla;
            
            r=[5 1 1000];
            [u, thetaV, gammaV, X_nfm,  U_nfm] = solve_nfm_over_time(o, X_target, T, thetaV_0,  o.cap*ones(1,T),r);
            %    theta= bsxfun(@rdivide, thetaV, sum(thetaV,1));
            
            h=figure;
            theta_plot = permute(thetaV(:,:,1:T),[1 3 2])' ;
            plot(1:T,theta_plot(:,1) ,'-bs',...
                1:T, theta_plot(:,2)  , '-ro',...
                1:T,  theta_plot(:,3)   , '--gd' ,...
                1:T,  theta_plot(:,4)   , '-.kd'); %   'LineWidth',2
            hleg1 = legend('host1,serv1','host2,serv1' ,'host3,serv1','host4,serv1');
            title('Transition of Service Placements');
            xlabel('Time'); ylabel('Portion of host multiplicity');
            UtilityLib.print_figure(h,9,7,sprintf('figure/sla_change_serv_portion' ));
            
            %                 h=figure;
            %                plot(1:T,U_nfm(1,:) ,'-bs',...
            %                     1:T, U_nfm(2,:)  , '-ro',...
            %                     1:T,  U_nfm(3,:)   , '--gd' ,...
            %                     1:T,  U_nfm(4,:)   , '-.kd');
            %                   hleg1 = legend('host1','host2' ,'host3','host4');
            %                  title('Host Utilizations');
            %                  xlabel('Time'); ylabel('Utilization');
            %                  UtilityLib.print_figure(h,9,7,sprintf('figure/sla_change_response_util' ));
            
            r=(N./X_nfm)-Z;
            h=figure;
            plot(1:T,r,'-bs');
            hleg1 = legend('R_{c1}');
            title('Class Response Time');
            xlabel('Time'); ylabel('Response Time');
            UtilityLib.print_figure(h,9,7,sprintf('figure/sla_change_class_response_time' ));
        end
        
        function test_change_sla_complex(o)
            T=100;
            o.setup_hosts();
            o.setup_services();
            r=[5 1 100];
            % o.c=rand(o.num_host,1)+10;
            N=[160
                100];
            Z=[1  1]';
            RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',sum(10)));
            c=rand(o.num_host,o.num_service,T);
            o.c = c(:,:,1);
            
            % before the sla change
            RT_sla0=[0.46  1.146 ]';
            x0=calculate_throuput_objective(o, N, Z, RT_sla0, 1);
            [theta,thetaV]=o.solve_nfm(x0,r);
            thetaV_0 =thetaV; % zeros(o.num_host , o.num_service);
            %print_lat_mat(o, thetaV)
            
            % after the sla change
            RT_sla1=[0.146  1.146 ]';
            x1=calculate_throuput_objective(o, N,  Z, RT_sla1, 1);
            [theta,thetaV]=o.solve_nfm(x1,r);
            [sum(thetaV_0,2) sum( thetaV,2) ]
            [sum(thetaV_0,1)' sum( thetaV,1)' ]
            % print_lat_mat(o, thetaV)
            
            
            % over time  
            N_t=N*ones(1,T);
            X_sla = [calculate_throuput_objective(o, N_t , Z, RT_sla0, T)...
                            calculate_throuput_objective(o, N_t , Z, RT_sla1, T)]; 
            X_target = X_sla;
            T=T*2;   
            o.c = repmat(c(:,:,1), [1,1,T]);   % o.c = c;
             [u, thetaV, gammaV, X_nfm,  U_nfm] = solve_nfm_over_time(o...
                 ,X_target, T, thetaV_0, o.cap*ones(1,T),r,'square');   
             
%             nsteps=T; 
%             T=7; 
%              N1=[N_t ones(1,T)*N_t(:,nsteps)];  
%             X_sla = calculate_throuput_objective(o, N1, Z, RT_sla, T+nsteps);
%             [u_all,thetaV_all, X_nfm_all, U_nfm_all,theta_all,c_infr, c_sla, c_trsh]=...
%                 iterate_simple_mpc_steps(o,X_sla,nsteps,T,r, thetaV_0,'square');
%             
            h=draw_thetaV2(o,thetaV,T,0.30);
            UtilityLib.print_figure(h,9,7,sprintf('figure/sla_change_theta' ));
            
            h=draw_hosts(o,U_nfm,T);
            UtilityLib.print_figure(h,9,7,sprintf('figure/sla_change_host_portion' ));
            
            
            h=draw_services(o,thetaV,T);
            UtilityLib.print_figure(h,9,7,sprintf('figure/sla_change_serv_portion' ));
            
            h=draw_classes_rt(o,X_nfm,X_target,N,N*ones(1,T),Z,T);
 %        h=draw_classes_rt(o,X_nfm_all,X_target,N,N*ones(1,T),Z,T);   %mpc 
            UtilityLib.print_figure(h,9,7,sprintf('figure/sla_change_rt' ));
            
        end
        
        
    end
    
end


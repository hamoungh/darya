classdef testmpc1_test_adding_service
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
             
        function test_adding_service(o)
            % o.setup_hosts();
            o.num_host=4;
            o.cap=[8 8 8 8]';
            o.speed_ratio=[1 1 1 1]';
            RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',sum(10)));
            o.cost=round(rand(o.num_host,1)*10)+1;
            
            % before change, steady state , initial state
            % after change
            o.num_classes = 2;
            o.num_service=2;
            call=zeros(o.num_service+o.num_classes,o.num_service+o.num_classes);
            call=[0 0 1 0
                0 0  0 1
                0 0 0  0
                0 0 0  0];
            o.c=rand(o.num_host,o.num_service)+10;
            o.c= [10.49851
                10.2248
                10.1981
                10.7605];
            o.call_graph = call;
            Y=calculate_contact_level(o,call);
            D =  [0.04  0.04];
            o.d= Y .* (ones(o.num_classes,1) * D);
            T=15;
            r=[5 1 1000];
            % before the service arrives
            N=[200 0]';
            Z=[1  1]';
            RT_sla=[0.146  1.146 ]';
            x0=calculate_throuput_objective(o, N(:,1), Z, RT_sla, 1);
            [theta,thetaV]=o.solve_nfm(x0,r);
            thetaV_0 =thetaV; % zeros(o.num_host , o.num_service);
            thetaV
            %print_lat_mat(o, thetaV)
            
            % after the service arrives
            N(2,1)=100; % N=[200 100]';
            x0=calculate_throuput_objective(o, N(:,1), Z, RT_sla, 1);
            [theta,thetaV]=o.solve_nfm(x0,r);
            thetaV
            % print_lat_mat(o, thetaV)
            
            % over time
            N_t=N*ones(1,T);
            X_sla = calculate_throuput_objective(o, N_t , Z, RT_sla, T);
            X_target = X_sla;
            
            %                 X_target = [X_sla
            %                                     zeros(1,T)];
            % % [sla violation cost,  infrastructure cost, change of
            % control cost]
            
            [u, thetaV, gammaV, X_nfm,  U_nfm] = solve_nfm_over_time(o...
                ,X_target, T, thetaV_0, o.cap*ones(1,T),r,'square');
            %    theta= bsxfun(@rdivide, thetaV, sum(thetaV,1));
            
            h=figure;
            theta_plot = permute(reshape(thetaV(2:3,:,1:T),[4,1,T]),[3,1,2]);
            plot(1:T,theta_plot(:,1) ,'-bs',...
                1:T, theta_plot(:,2)  , '-ro',...
                1:T,  theta_plot(:,3)   , '--gd' ,...
                1:T,  theta_plot(:,4)   , '-.kd'); %   'LineWidth',2
            hleg1 = legend('host2,serv1','host3,serv1' ,'host2,serv2','host3,serv2');
            title('Transition of Service Placements');
            xlabel('Time'); ylabel('Portion of host multiplicity');
            UtilityLib.print_figure(h,9,7,sprintf('figure/adding_class_response_serv_portion' ));
            
            h=figure;
            plot(1:T,U_nfm(1,:) ,'-bs',...
                1:T, U_nfm(2,:)  , '-ro',...
                1:T,  U_nfm(3,:)   , '--gd' ,...
                1:T,  U_nfm(4,:)   , '-.kd');
            hleg1 = legend('host1','host2' ,'host3','host4');
            title('Host Utilizations');
            xlabel('Time'); ylabel('Utilization');
            UtilityLib.print_figure(h,9,7,sprintf('figure/adding_class_response_util' ));
            
            h=figure;
            plot(1:T,X_nfm(1,:),'-bs', 1:T,X_target(1,:),'--b',  1:T, X_nfm(2,:) ,'-ro', 1:T,X_target(2,:),'--r' );
            hleg1 = legend('X_{c1}', 'X_{c1}^{SLA}', 'X_{c2}', 'X_{c2}^{SLA}');
            title('Class Througput');
            xlabel('Time'); ylabel('Througput');
            UtilityLib.print_figure(h,9,7,sprintf('figure/adding_class_throughput' ));
        end
        
        
        function test_adding_service_complex (o)
            T=100;
            o.setup_hosts();
            o.setup_services();
            r=[5 1 100];
            % o.c=rand(o.num_host,1)+10;
            N=[160
                100];
            Z=[1  1]';
            RT_sla=[0.146  1.146 ]';
            x0=calculate_throuput_objective(o, N(:,1), Z, RT_sla, 1);
            
            %                 c=...
            %                     [10.1981
            %                    10.7605
            %                    10.1691
            %                    10.0883
            %                    10.6854
            %                    10.9534]*ones(1,o.num_service);
            
            RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',sum(10)));
            c=rand(o.num_host,o.num_service,T);
            o.c = c(:,:,1);
            % before the service arrives
            N1=N; N1(2,1)=0;
            x0=calculate_throuput_objective(o, N1, Z, RT_sla, 1);
            [theta,thetaV]=o.solve_nfm(x0,r);
            thetaV_0 =thetaV; % zeros(o.num_host , o.num_service);
            %print_lat_mat(o, thetaV)
            
            % after the service arrives
            x0=calculate_throuput_objective(o, N,  Z, RT_sla, 1);
            [theta,thetaV]=o.solve_nfm(x0,r);
            [sum(thetaV_0,2) sum( thetaV,2) ]
            [sum(thetaV_0,1)' sum( thetaV,1)' ]
            % print_lat_mat(o, thetaV)
            
            % over time
            N_t=N*ones(1,T);
            o.c = repmat(c(:,:,1), [1,1,T]);   % o.c = c;
            X_sla = calculate_throuput_objective(o, N_t , Z, RT_sla, T);
            X_target = X_sla;
            [u, thetaV, gammaV, X_nfm,  U_nfm] = solve_nfm_over_time(o...
                ,X_target, T, thetaV_0, o.cap*ones(1,T),r,'square');
            
            h=draw_hosts(o,U_nfm,T);
            UtilityLib.print_figure(h,9,7,sprintf('figure/adding_class_host_portion' ));
            
            h=draw_services(o,thetaV,T);
            UtilityLib.print_figure(h,9,7,sprintf('figure/adding_class_serv_portion' ));
            
            h=draw_classes_throughput(o,X_nfm,X_target,T);
            UtilityLib.print_figure(h,9,7,sprintf('figure/adding_class_throughput' ));
            
            win=30;
            h=draw_thetaV(o,thetaV,win,T);
            UtilityLib.print_figure(h,9,7,sprintf('figure/adding_class_theta' ));
            % status = close(h)
            
        end
        
        function test_adding_unmanaged_service(o)
            % o.setup_hosts();
            o.num_host=4;
            o.cap=[8 8 8 8]';
            o.speed_ratio=[1 1 1 1]';
            RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',sum(10)));
            o.cost=round(rand(o.num_host,1)*10)+1;
            
            % before change, steady state , initial state
            % after change
            o.num_classes = 2;
            o.num_service=2;
            call=zeros(o.num_service+o.num_classes,o.num_service+o.num_classes);
            call=[0 0 1 0
                0 0  0 1
                0 0 0  0
                0 0 0  0];
            o.c=rand(o.num_host,o.num_service)+10;
            o.c= [10.49851
                10.2248
                10.1981
                10.7605];
            o.c=ones(4,1);
            o.call_graph = call;
            Y=calculate_contact_level(o,call);
            D =  [0.04  0.04];
            o.d= Y .* (ones(o.num_classes,1) * D);
            T=15;
            
            % initial steady state: before the service arrives
            N=[0
                400];
            Z=[1  1]';
            RT_sla=[0.546  1.146 ]';
            x0=calculate_throuput_objective(o, N(:,1), Z, RT_sla, 1);
            X_target=x0;
            d=o.d';
            r1 = 1000*ones(1,o.num_classes);
            cvx_begin quiet
            variables  thetaV(o.num_host,o.num_service)...
                gammaV(o.num_service,o.num_classes) ...
                X_c(o.num_classes,1)
            thetaV>=0;
            thetaV(:,1).*[1 0 0 0]' == thetaV(:,1);
            sum(thetaV,2)<=o.cap .* o.speed_ratio
            sum(thetaV,1)' == sum(gammaV ,2)
            gammaV == d .* (ones(o.num_service,1)*X_c')
            X_target <= X_c;
            minimize sum( sum( thetaV,2))
            cvx_end
            thetaV_0 =thetaV; gammaV_0=gammaV;
            thetaV
            % print_lat_mat(o, thetaV)
            
            % final steady state
            N=[300
                400];
            Z=[1  1]';
            RT_sla=[0.546  1.146 ]';
            x0=calculate_throuput_objective(o, N(:,1), Z, RT_sla, 1);
            X_target=x0;
            d=o.d';
            r1 = 1000*ones(1,o.num_classes);
            cvx_begin quiet
            variables  thetaV(o.num_host,o.num_service)...
                gammaV(o.num_service,o.num_classes) ...
                X_c(o.num_classes,1)
            thetaV>=0;
            thetaV(:,1).*[1 0 0 0]' == thetaV(:,1);
            sum(thetaV,2)<=o.cap .* o.speed_ratio
            sum(thetaV,1)' == sum(gammaV ,2)
            gammaV == d .* (ones(o.num_service,1)*X_c')
            X_target <= X_c;
            minimize sum( sum( thetaV,2))
            cvx_end
            thetaV_end = thetaV;
            thetaV
            % print_lat_mat(o, thetaV)
            
            % transition over time
            N=[ones(1,T)*300
                ones(1,T)*400];
            X_sla = calculate_throuput_objective(o, N, Z, RT_sla, T);
            f_sla = X_sla;
            r=[5 1 700];
            %   [u, thetaV, gammaV, X_nfm,  U_nfm] = solve_nfm_over_time(o, X_target, T, thetaV_0, o.cap*ones(1,T),r);
            cap=o.cap*ones(1,T);
            % function [u, thetaV,gammaV, f_c, U_host_, theta ]=solve_nfm_over_time(o, f_sla, T, thetaV_0, cap, r)
            d=o.d';
            cvx_begin  quiet
            cvx_precision best
            variables  thetaV(o.num_host , o.num_service, T+1) ... % state variable
                u(o.num_host , o.num_service, T) ...    % input
                gammaV(o.num_service , o.num_classes,T) ... % feed forward input variable, directly controlable
                f_c(o.num_classes,T)...  %  output variable in MPC  based on input
                U_host(o.num_host,T) ... % output, utilization
                consumption_cost(o.num_host,T) ... % cost of output
                consumption_cost_tot(1,T)       % cost of output
            thetaV>=0
            
            thetaV(:,1,:).*repmat([1 0 0 0]',[1 1 T+1]) == thetaV(:,1,:);
            thetaV(:,:,1) ==  thetaV_0
            thetaV(:,:,2:T+1) ==  thetaV(:,:,1:T) + u(:,:,1:T)
            U_host -  cap  <= 0 ;
            
            permute(sum(thetaV(:,:,1:T),1),[2,3,1])  == permute(sum(gammaV, 2),[1,3,2])
            gammaV(:,:,1:T) == ...
                repmat(d,[1,1,(T)]) .* reshape(...
                (ones(o.num_service,1)* reshape(f_c(:,1:T),[1,o.num_classes*(T)])), ...
                [o.num_service,o.num_classes, (T)])
            U_host == permute(sum(thetaV(:,:, 1:T),2),[1,3,2]) ;
            consumption_cost_tot ==  o.c'* U_host
            minimize(... % time
                r(1) *  sum(sum(pos( f_sla-f_c(:,1:T)),1),2)   +  ...  % sla violation cost
                r(2) *  sum(consumption_cost_tot(:,1:T)) + ...      % infrastructure cost
                r(3) *  sum(permute(sum(square(u),1),[3,2,1])*[0 1]')  )
            %  r(3) *  sum(sum(sum(square(u),2),1),3) )       % u=o.num_host , o.num_service, T
            cvx_end
            U_host_=U_host;
            theta= bsxfun(@rdivide, thetaV, sum(thetaV,1));
            h=figure;
            theta_plot = permute(thetaV(1,1,1:T),[3,1,2]);
            plot(1:T,permute(thetaV(1,1,1:T),[3,1,2]) ,'-bs',...
                1:T, permute(thetaV(1,2,1:T),[3,1,2]) , '-ro',...
                1:T, permute(thetaV(2,2,1:T),[3,1,2])   , '--gd' ); %   'LineWidth',2
            hleg1 = legend('host1,serv1','host1,serv2' ,sprintf('host2,serv2\n host3,serv2\n host4,serv2'));
            title('Transition of Service Placements');
            xlabel('Time'); ylabel('Portion of host multiplicity');
            UtilityLib.print_figure(h,9,7,sprintf('figure/migrate_on_external_class_vtheta' ));
            
            h=figure;
            plot(1:T,U_host(1,:) ,'-bs',...
                1:T, U_host(2,:)  , '-ro',...
                1:T,  U_host(3,:)   , '--gd' ,...
                1:T,  U_host(4,:)   , '-.kd');
            hleg1 = legend('host1','host2' ,'host3','host4');
            title('Host Utilizations');
            xlabel('Time'); ylabel('Utilization');
            UtilityLib.print_figure(h,9,7,sprintf('figure/migrate_on_external_class_util' ));
            
            %                  h=figure;
            %                  X_nfm=f_c;
            %                  plot(1:T,X_nfm(1,:),'-bs', 1:T,X_target(1,:),'--b',  1:T, X_nfm(2,:) ,'-ro', 1:T,X_target(2,:),'--r' );
            %                   hleg1 = legend('X_{c1}', 'X_{c1}^{SLA}', 'X_{c2}', 'X_{c2}^{SLA}');
            %                  title('Class Througput');
            %                  xlabel('Time'); ylabel('Througput');
            %                  UtilityLib.print_figure(h,9,7,sprintf('figure/migrate_on_external_class_throughput' ));
            
        end
        
        function test_adding_unmanaged_service_complex (o)
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
            RT_sla=[0.46  1.146 ]';
            
            
            % initial steady state: before the service arrives
            N1=N; N1(1,1)=0;
            N=[0
                400];
            Z=[1  1]';
            x0=calculate_throuput_objective(o, N(:,1), Z, RT_sla, 1);
            X_target=x0;
            d=o.d';
            cvx_begin quiet
            variables  thetaV(o.num_host,o.num_service)...
                gammaV(o.num_service,o.num_classes) ...
                X_c(o.num_classes,1)
            thetaV>=0;
            thetaV(:,1).*[1 0 0 0]' == thetaV(:,1);
            sum(thetaV,2)<=o.cap .* o.speed_ratio
            sum(thetaV,1)' == sum(gammaV ,2)
            gammaV == d .* (ones(o.num_service,1)*X_c')
            X_target <= X_c;
            minimize sum( sum(o.c .* thetaV,2),1);  % sum( sum( thetaV,2))
            cvx_end
            thetaV_0 =thetaV; gammaV_0=gammaV;
            thetaV
            % print_lat_mat(o, thetaV)
            
            % final steady state
            N=[300
                400];
            Z=[1  1]';
            RT_sla=[0.546  1.146 ]';
            x0=calculate_throuput_objective(o, N(:,1), Z, RT_sla, 1);
            X_target=x0;
            d=o.d';
            r1 = 1000*ones(1,o.num_classes);
            cvx_begin quiet
            variables  thetaV(o.num_host,o.num_service)...
                gammaV(o.num_service,o.num_classes) ...
                X_c(o.num_classes,1)
            thetaV>=0;
            thetaV(:,1).*[1 0 0 0]' == thetaV(:,1);
            sum(thetaV,2)<=o.cap .* o.speed_ratio
            sum(thetaV,1)' == sum(gammaV ,2)
            gammaV == d .* (ones(o.num_service,1)*X_c')
            X_target <= X_c;
            minimize sum( sum(o.c .* thetaV,2),1);   % sum( sum( thetaV,2))
            cvx_end
            thetaV_end = thetaV;
            thetaV
            % print_lat_mat(o, thetaV)
            
            % transition over time
            N=[ones(1,T)*300
                ones(1,T)*400];
            X_sla = calculate_throuput_objective(o, N, Z, RT_sla, T);
            f_sla = X_sla;
            r=[5 1 700];
            %   [u, thetaV, gammaV, X_nfm,  U_nfm] = solve_nfm_over_time(o, X_target, T, thetaV_0, o.cap*ones(1,T),r);
            cap=o.cap*ones(1,T);
            % function [u, thetaV,gammaV, f_c, U_host_, theta ]=solve_nfm_over_time(o, f_sla, T, thetaV_0, cap, r)
            d=o.d';
            cvx_begin  quiet
            cvx_precision best
            variables  thetaV(o.num_host , o.num_service, T+1) ... % state variable
                u(o.num_host , o.num_service, T) ...    % input
                gammaV(o.num_service , o.num_classes,T) ... % feed forward input variable, directly controlable
                f_c(o.num_classes,T)...  %  output variable in MPC  based on input
                U_host(o.num_host,T) ... % output, utilization
                consumption_cost(o.num_host,T) ... % cost of output
                consumption_cost_tot(1,T)       % cost of output
            thetaV>=0
            
            thetaV(:,1,:).*repmat([1 0 0 0]',[1 1 T+1]) == thetaV(:,1,:);
            thetaV(:,:,1) ==  thetaV_0
            thetaV(:,:,2:T+1) ==  thetaV(:,:,1:T) + u(:,:,1:T)
            U_host -  cap  <= 0 ;
            
            permute(sum(thetaV(:,:,1:T),1),[2,3,1])  == permute(sum(gammaV, 2),[1,3,2])
            gammaV(:,:,1:T) == ...
                repmat(d,[1,1,(T)]) .* reshape(...
                (ones(o.num_service,1)* reshape(f_c(:,1:T),[1,o.num_classes*(T)])), ...
                [o.num_service,o.num_classes, (T)])
            U_host == permute(sum(thetaV(:,:, 1:T),2),[1,3,2]) ;
            consumption_cost_tot ==  o.c'* U_host
            minimize(... % time
                r(1) *  sum(sum(pos( f_sla-f_c(:,1:T)),1),2)   +  ...  % sla violation cost
                r(2) *  sum(consumption_cost_tot(:,1:T)) + ...      % infrastructure cost
                r(3) *  sum(permute(sum(square(u),1),[3,2,1])*[0 1]')  )
            %  r(3) *  sum(sum(sum(square(u),2),1),3) )       % u=o.num_host , o.num_service, T
            cvx_end
            U_host_=U_host;
            theta= bsxfun(@rdivide, thetaV, sum(thetaV,1));
            h=figure;
            theta_plot = permute(thetaV(1,1,1:T),[3,1,2]);
            plot(1:T,permute(thetaV(1,1,1:T),[3,1,2]) ,'-bs',...
                1:T, permute(thetaV(1,2,1:T),[3,1,2]) , '-ro',...
                1:T, permute(thetaV(2,2,1:T),[3,1,2])   , '--gd' ); %   'LineWidth',2
            hleg1 = legend('host1,serv1','host1,serv2' ,sprintf('host2,serv2\n host3,serv2\n host4,serv2'));
            title('Transition of Service Placements');
            xlabel('Time'); ylabel('Portion of host multiplicity');
            UtilityLib.print_figure(h,9,7,sprintf('figure/migrate_on_external_class_vtheta' ));
            
            h=figure;
            plot(1:T,U_host(1,:) ,'-bs',...
                1:T, U_host(2,:)  , '-ro',...
                1:T,  U_host(3,:)   , '--gd' ,...
                1:T,  U_host(4,:)   , '-.kd');
            hleg1 = legend('host1','host2' ,'host3','host4');
            title('Host Utilizations');
            xlabel('Time'); ylabel('Utilization');
            UtilityLib.print_figure(h,9,7,sprintf('figure/migrate_on_external_class_util' ));
            
            %                  h=figure;
            %                  X_nfm=f_c;
            %                  plot(1:T,X_nfm(1,:),'-bs', 1:T,X_target(1,:),'--b',  1:T, X_nfm(2,:) ,'-ro', 1:T,X_target(2,:),'--r' );
            %                   hleg1 = legend('X_{c1}', 'X_{c1}^{SLA}', 'X_{c2}', 'X_{c2}^{SLA}');
            %                  title('Class Througput');
            %                  xlabel('Time'); ylabel('Througput');
            %                  UtilityLib.print_figure(h,9,7,sprintf('figure/migrate_on_external_class_throughput' ));
            
        end
 
    end
    
end


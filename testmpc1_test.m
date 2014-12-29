classdef testmpc1_test
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
            % ----------------------------------------------------------------------------
        %  Methods  for testing different cases
        % ----------------------------------------------------------------------------
          
        % note: use   r(3) *  sum(sum_square(sum(u,2),1),3) )  in the
        % cost function
        function test_change_cost_coef(o)
            T=100;
            o.setup_hosts();
            o.setup_services();
            r=[5 1 10];
            % o.c=rand(o.num_host,1)+10;
            N=[160
                100];
            Z=[1  1]';
            RT_sla=[0.146  1.146 ]';
            x0=calculate_throuput_objective(o, N(:,1), Z, RT_sla, 1);
            x0
            % before
            %                  o.c=[10.1981
            %                    10.7605
            %                    10.1691
            %                    10.0883
            %                    10.6854
            %                    10.9534]*ones(1,o.num_service);
            RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',sum(10)));
            o.c=rand(o.num_host,o.num_service);
            [B, IX] = sort(o.c);
            [theta,thetaV]=o.solve_nfm(x0,r);
            thetaV_0 =thetaV;
            thetaV
            spy(thetaV,'bs',10); hold on;
            
            util_x0=sum(thetaV,2)
            % print_lat_mat(o, thetaV)
            
            % after the service arrives
            %                 o.c=[10.9039
            %                        10.5122
            %                        10.8126
            %                        10.9125
            %                        10.7218
            %                        10.2919]*ones(1,o.num_service);
            RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',sum(1110)));
            o.c=rand(o.num_host,o.num_service);
            %  print_lat_mat(o, thetaV)
            [B, IX] = sort(o.c);
            [theta,thetaV]=o.solve_nfm(x0,r);
            thetaV
            spy(thetaV,'rs',10)
            util_xn=sum(thetaV,2)
            
            % over time
            N_t=N*ones(1,T);
            X_sla = calculate_throuput_objective(o, N_t, Z, RT_sla, T);
            X_target = X_sla;
            
            
            [u, thetaV, gammaV, X_nfm,  U_nfm] = solve_nfm_over_time(o, X_target, T, thetaV_0, o.cap*ones(1,T),r, 'square');
            thetaV(:,:,end-1)
            %  u(:,:,1)
            [U_nfm(:,1) U_nfm(:,end)]
            
            h=draw_hosts(o,U_nfm,T)
            UtilityLib.print_figure(h,9,7,sprintf('figure/test_change_cost_coef_util' ));
            
            trsh = 0.15;
            h=draw_thetaV2(o,thetaV,T,trsh);
            UtilityLib.print_figure(h,9,7,sprintf('figure/test_change_cost_coef_theta' ));
            
            win=30;
            h=draw_thetaV(o,thetaV,win,T);
            UtilityLib.print_figure(h,9,7,sprintf('figure/test_change_cost_coef_theta_matrix' ));
        end
  
        function test_nfm_static(o)
            % setup_hosts(o)
            o.num_host=4;
            o.cap=[8 8 8 8]';
            o.speed_ratio=[1 1 1 1]';
            RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',sum(10)));
            o.cost=round(rand(o.num_host,1)*10)+1;
            
            % setup_services
            o.num_classes = 1;
            o.num_service=1;
            call=zeros(o.num_service+o.num_classes,o.num_service+o.num_classes);
            call=[1 1
                1 1];
            o.c=rand(o.num_host,o.num_service)+10;
            o.c= [10.49851
                10.2248
                10.1981
                10.7605];
            o.c;
            
            o.call_graph = call;
            Y=calculate_contact_level(o,call);
            
            D =  [0.04 ];
            o.d= Y .* (ones(o.num_classes,1) * D);
            
            % test
            N=[250]';
            Z=[1]';
            RT_sla=[0.1]';
            
            
            % initialize the loop
            X_sla =  N./(RT_sla + Z) ;  % here each element of RT is \sum_k RT_c,k
            X_target = X_sla;
            test=true;
            MAX_ITER= 20;
            iter = 1;
            disp(sprintf('X_sla=%4.3f ', X_sla ));
            while test & iter<MAX_ITER
                disp(sprintf('---------------------- iteration %d----------------- ', iter));
                disp(sprintf('X_target=%4.3f',  X_target ));
                [theta]= solve_nfm(o, X_target);   % thetaV(o.num_host,o.num_service)
                
                disp('invocations');    %  disp(sprintf('invocations are  %20.3f', theta ));
                disp(theta');
                
                % achieved throughput, which might be even different from sla
                % X_nfm = gammaMu ./ o.d;
                
                X_target_old = X_target;
                % d_csh= permute(  repmat(o.d, [ 1,1 , o.num_host] ) , [2,1,3])  .* ...
                %             permute( repmat(theta, [1 , 1, o.num_classes]), [3,2,1]  );
                
                %              [X_q,R_q] = solve_lqm2(o,d_csh,N,Z);
                %            disp( sprintf('octave says %20.3f %20.3f', X_q,R_q ));
                
                [X_q,R_q] =  solve_lqm(o, theta, N, Z);
                disp(sprintf('X_q=%4.3f \n R_q=%4.3f\n', X_q,R_q ));
                X_target = X_target + (X_sla - X_q);
                
                if ((X_target-X_target_old) < 0.001*X_target)  break; end;
                iter=iter+1;
            end
            
            % solve_lqm(o,d,N,Z,)
        end
        
        % depricated
        %  assumption of controllability
        % successive approximation method, which makes multiple calls to the underlying solver
        function test_nfm_mpc_1class(o)
            % o.setup_hosts();
            o.num_host=4;
            o.cap=[8 8 8 8]';
            o.speed_ratio=[1 1 1 1]';
            RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',sum(10)));
            o.cost=round(rand(o.num_host,1)*10)+1;
            
            % setup_services
            o.num_classes = 1;
            o.num_service=1;
            call=zeros(o.num_service+o.num_classes,o.num_service+o.num_classes);
            call=[1 1
                1 1];
            o.c=rand(o.num_host,o.num_service)+10;
            o.c= [10.49851
                10.2248
                10.1981
                10.7605];
            
            o.call_graph = call;
            Y=calculate_contact_level(o,call);
            
            D =  [0.04 ];
            o.d= Y .* (ones(o.num_classes,1) * D);
            
            T=3;
            %  o.setup_services();
            N=[250 100 200];
            Z=[1]';
            RT_sla=[0.146]';
            %   o.r1 = 1000*ones(1,o.num_classes);
            
            thetaV_0 = zeros(o.num_host , o.num_service);
            
            X_sla = calculate_throuput_objective(o, N, Z, RT_sla, T);
            
            
            X_target = X_sla;
            test=true;
            MAX_ITER= 20;
            iter = 1;
            thetaV=[];
            disp(sprintf('X_sla=%4.3f ', X_sla ));
            while test & iter<MAX_ITER
                disp(sprintf('---------------------- iteration %d----------------- ', iter));
                disp('X_target' );
                X_target
                %   [theta]= solve_nfm(o, X_target);   % thetaV(o.num_host,o.num_service)
                [u, thetaV, gammaV, X_nfm,  U_nfm] = solve_nfm_over_time(o, X_target, T, thetaV_0, o.cap*ones(1,T+1));
                
                %  disp('invocations');
                disp( thetaV(:,:, 1:T+1) );
                
                X_target_old = X_target;
                % d_csh= permute(  repmat(o.d, [ 1,1 , o.num_host] ) , [2,1,3])  .* ...
                %             permute( repmat(theta, [1 , 1, o.num_classes]), [3,2,1]  );
                
                % solve the model for each time step
                X_q = []; R_q=[];
                for t=1:T
                    [x_q, r_q] =  solve_lqm(o, thetaV(:,:,t+1)  , N(:,t),  Z);
                    X_q = [X_q x_q];
                    R_q = [R_q r_q];
                end
                
                disp('X_q');  disp(X_q);
                disp('R_q');  disp(R_q);
                X_target = X_target + (X_sla - X_q);
                disp(sprintf('X_sla - X_q=%4.3f\n', X_sla - X_q ));
                if ( (X_target-X_target_old) < 0.001*X_target)  break; end;
                iter=iter+1;
            end
        end
        
        %  assumption of controllability
        % successive approximation method, which makes multiple calls to the underlying solver
        function test_nfm_mpc_2class(o)
            % o.setup_hosts();
            o.num_host=4;
            o.cap=[8 8 8 8]';
            o.speed_ratio=[1 1 1 1]';
            RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',sum(10)));
            o.cost=round(rand(o.num_host,1)*10)+1;
            
            % setup_services
            o.num_classes = 2;
            o.num_service=1;
            % setting  the call multiplicity from all centers including
            % delay+queuing center to one another
            % here both classes call the service
            call=zeros(o.num_service+o.num_classes,o.num_service+o.num_classes);
            call=[0 0 1
                0 0  1
                0 0 0];
            o.call_graph = call;
            Y=calculate_contact_level(o,call);
            % D is just services
            D =  [0.04 ];
            % d is class*service
            o.d= Y .* (ones(o.num_classes,1) * D);
            
            o.c=rand(o.num_host,o.num_service)+10;
            o.c= [10.49851
                10.2248
                10.1981
                10.7605];
            
            T=3;
            %  o.setup_services();
            N=[250 100 200
                100   250 50];
            Z=[1 1]';
            RT_sla=[0.146  0.267]';
            objective='abs'; 
            r=[50 1 1];  
            thetaV_0=zeros(o.num_host , o.num_service);  
            X_sla = calculate_throuput_objective(o, N, Z, RT_sla, T);
                
            thetaV= solve_mpc_contention...%(o, T, N, Z, RT_sla )
                (o,N, Z, X_sla, T, thetaV_0,o.cap*ones(1,T), r, objective);
        end
    
    
        
        % [optcost_tr, costmpc, optcost] = testmpc1().test_dynamic_workload(50, 5, 400, 'square')
        function [optcost_tr, mpccost, optcost] = test_dynamic_workload(o, nsteps, T,r_trsh_cost, objective)
            o.num_host=4;
            o.num_classes = 1;
            o.num_service=1;
            o.cap=[8 8 8 8]';
            o.speed_ratio=[1 1 1 1]';
            RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',10));
            % o.cost=round(rand(o.num_host,1)*10)+1;
            o.c=rand(o.num_host,1)*ones(1, o.num_service)+10;
            
            call=zeros(o.num_service+o.num_classes,o.num_service+o.num_classes);
            call=[1 1
                1 1];
            o.call_graph = call;
            Y=calculate_contact_level(o,call);
            D =  [0.04 ];
            o.d= Y .* (ones(o.num_classes,1) * D);
            Z=[1]';
            RT_sla=[0.146]';
             r=[50 1 r_trsh_cost];
             
            %             % mpc d
            % nsteps = 10;
            %  T=5;
            thetaV_0 =   [0   8.0000   0.0360  0]';
           
            N=200+[1:nsteps]+50*sin([1:5:5*(nsteps)]/30); 
            N1=[N ones(1,T)*N(:,nsteps)];
            h=figure; 
             plot(1:nsteps, N ,'-r' );  
            hleg1 = legend('N');
            title('Workload');
            xlabel('Time step'); ylabel('Number of Users (N)');  
        %    axis([0 nsteps min(N_orig)-20 max(N_orig)+20])
            UtilityLib.print_figure(h,9,7,sprintf('figure/dynamic_workload'));
            % plot(N); hold on;
            X_sla = calculate_throuput_objective(o, N1, Z, RT_sla, T+nsteps);
           
            [u_all,thetaV_all, X_nfm_all, U_nfm_all,theta_all,c_infr, c_sla, c_trsh]=...
                iterate_simple_mpc_steps(o,X_sla,nsteps,T,r, thetaV_0,objective);
            disp(sprintf('\n c_sla=%5.4f,  c_infr=%5.4f, c_trsh=%5.4f', c_sla,  c_infr, c_trsh));
            mpccost =  r(1)*c_sla + r(2) * c_infr+ r(3)*c_trsh
            permute(sum(abs(u_all(:,:,:))), [1 3 2]) ;
            
            % this one igores the cost of trashing totally,
            % so it is off
            T_orig=T;
            T=nsteps;
            % N=200+[1:T]+50*sin([1:5:5*T]/10);
            X_sla = calculate_throuput_objective(o, N, Z, RT_sla, T);
            % thetaV_0 = zeros(o.num_host , o.num_service);
            thetaV_0 =   [0   8.0000   0.0360  0]';
            X_target = X_sla;
            r_trsh=[50 1 0.01];
            [u_opt_tr, thetaV, gammaV, X_nfm,  U_nfm_tr, theta] =...
                solve_nfm_over_time(o, X_target, T, thetaV_0, o.cap*ones(1,T), r_trsh, objective);
            c_infr_opt_tr=sum(o.c' * U_nfm_tr);
            % we do not evaluate for the steps we do not control
            % [X_lqm ~]=eveluate_lqm(o, theta , N,  Z);
            c_sla_opt_tr= sum(pos( X_sla(:,2:T)-X_nfm(:,2:T)));
            c_trsh_opt_tr= sum(sum(sum(abs(u_opt_tr(:,:,:)))));
            disp(sprintf('\n c_sla_opt_tr=%5.4f,  c_infr_opt_tr=%5.4f,  c_trsh_opt_tr=%5.4f,', ...
                c_sla_opt_tr,  c_infr_opt_tr,  c_trsh_opt_tr));
            optcost_tr =  r(1)*c_sla_opt_tr + r(2) * c_infr_opt_tr + r(3)*c_trsh_opt_tr
            permute(sum(abs(u_opt_tr(:,:,:))), [1 3 2]);
            
            % this is the optimal for the simplistic approach evaluated
            % with lqm
            [u_opt, thetaV, gammaV, X_nfm,  U_nfm_notr, theta] =...
                solve_nfm_over_time(o, X_target, T, thetaV_0, o.cap*ones(1,T), r, objective);
            permute(u_opt,[1 3 2]);
            c_infr_opt=sum(o.c' * U_nfm_notr);
            % we do not evaluate for the steps we do not control
            % c_sla_opt= sum(pos( X_sla(:,2:T)-X_nfm(:,2:T)));
            % [X_lqm ~]=eveluate_lqm(o, theta , N,  Z);
            c_sla_opt= sum(pos( X_sla(:,2:T)-X_nfm(:,2:T)));
            c_trsh_opt= sum(sum(sum(abs(u_opt(:,:,:)))));
            disp(sprintf('\n c_sla_opt=%5.4f, c_infr_opt=%5.4f, c_trsh_opt=%5.4f',...
                c_sla_opt,  c_infr_opt,  c_trsh_opt));
            optcost =  r(1)*c_sla_opt + r(2) * c_infr_opt + r(3)*c_trsh_opt
            
            % compare cost of infrastructure
            %  plot(1:T, o.c' *U_nfm_tr,'b-', 1:T, o.c' *U_nfm_notr ,'r--' );
            h=figure;
            plot(1:T, o.c' *U_nfm_tr,'-bs',...
                    1:T, o.c' *U_nfm_notr ,'-ro',... 
                    1:T, o.c' *U_nfm_all, '--gd'); 
            hleg1 = legend('trashing ignored','optimal' ,'mpc');
            title('Cost of infrastructure');
            xlabel('Time step'); ylabel('Infrastructure Cost');
            UtilityLib.print_figure(h,9,7,sprintf('figure/dynamic_infrastructure_cost_nsteps%d_T%d_rTRSH%d%s', nsteps, T_orig,r_trsh_cost,objective));
            % compare cost of trashing
            h=figure;
            plot(1:T, permute(sum(abs(u_opt_tr(:,:,:))), [1 3 2]),  '-bs',...
                1:T,  permute(sum(abs(u_opt(:,:,:))), [1 3 2]),'-ro',...
                1:T, permute(sum(abs(u_all(:,:,:))), [1 3 2])  , '--gd');
            hleg1 = legend('trashing ignored','optimal' ,'mpc');
            title('Cost of trashing');
            xlabel('Time step'); ylabel('Trashing Cost');
            UtilityLib.print_figure(h,9,7,sprintf('figure/dynamic_trashing_cost_nsteps%d_T%d_rTRSH%d%s', nsteps, T_orig,r_trsh_cost,objective));
        end
        
        % testmpc1().draw_cost_with_mpc_horizon( 10 , 40, 5, 7)
        function draw_cost_with_mpc_horizon(o,nsteps, r_trsh_cost, horizon_lb, horizon_ub)
            o.num_host=4;
            o.num_classes = 1;
            o.num_service=1;
            o.cap=[8 8 8 8]';
            o.speed_ratio=[1 1 1 1]';
            RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',10));
            % o.cost=round(rand(o.num_host,1)*10)+1;
            o.c=rand(o.num_host,1)*ones(1, o.num_service)+10;
            
            %                  o.c=rand(o.num_host,o.num_service)+10;
            %                  o.c= [10.49851
            %                         10.2248
            %                         10.1981
            %                         10.7605];
            call=zeros(o.num_service+o.num_classes,o.num_service+o.num_classes);
            call=[1 1
                1 1];
            o.call_graph = call;
            Y=calculate_contact_level(o,call);
            D =  [0.04 ];
            o.d= Y .* (ones(o.num_classes,1) * D);
            Z=[1]';
            RT_sla=[0.146]';
            
            thetaV_0 =   [0   8.0000   0.0360  0]';
            
            % plot(N); hold on;
            
            r=[50 1 r_trsh_cost];
            cost=[];
            N=200+[1:nsteps]+50*sin([1:5:5*(nsteps)]/10);
            for horizon=horizon_lb:horizon_ub
                T=horizon;
                N_=[N ones(1,T)*N(nsteps)];
                X_sla = calculate_throuput_objective(o, N_, Z, RT_sla, T+nsteps);
                
                [u_all, thetaV_all, X_nfm_all, U_nfm_all,theta_all,c_infr, c_sla, c_trsh]= iterate_simple_mpc_steps(o,X_sla,nsteps,T,r, thetaV_0);
                mpccost =  r(1)*c_sla + r(2) * c_infr+ r(3)*c_trsh;
                fprintf('%4.3f, %4.3f, %4.3f, %4.3f',c_sla,  c_infr, c_trsh, mpccost);
                cost=[cost; r(1)*c_sla, r(2) *  c_infr,  r(3)*c_trsh, mpccost];
                permute(sum(abs(u_all(:,:,:))), [1 3 2]) ;
            end
            
            h=figure;
            plot(horizon_lb:horizon_ub, cost(:,4) ,'-bs', ...
                horizon_lb:horizon_ub, cost(:,2) ,'-ro',...
                horizon_lb:horizon_ub, cost(:,3) , '--gd');
            hleg1 = legend('total cost',sprintf('%d*infrastructure cost',r(2)) , sprintf('%d*trashing cost', r(3) ));
            title('Cost over Lookahead Horizon');
            xlabel('Lookahead Horizon (N)'); ylabel('Cost');
            UtilityLib.print_figure(h,9,7,sprintf('figure/cost_over_lookahead_nsteps%d_rTrshCost%', nsteps,r_trsh_cost));
        end
        
        function test_draw_tradeoff_curve(o)
            o.num_host=4;
            o.num_classes = 1;
            o.num_service=1;
            o.cap=[8 8 8 8]';
            o.speed_ratio=[1 1 1 1]';
            RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',10));
            o.c=rand(o.num_host,1)*ones(1, o.num_service)+10;
            
            call=zeros(o.num_service+o.num_classes,o.num_service+o.num_classes);
            call=[1 1
                1 1];
            o.call_graph = call;
            Y=calculate_contact_level(o,call);
            D =  [0.04 ];
            o.d= Y .* (ones(o.num_classes,1) * D);
            Z=[1]';
            RT_sla=[0.146]';
            
            cost_infr=[]; cost_trsh=[];
            r_size=40;
            % r_=logspace( 0, 2, r_size )
            r_=linspace(1,100,r_size);
            r=[50 1 0];
            thetaV_0 =   [0   8.0000   0.0360  0]';
            T=100;
            N=200+[1:T]+50*sin([1:5:5*T]/10);
            X_sla = calculate_throuput_objective(o, N, Z, RT_sla, T);
            X_target=X_sla;
            for i=1:size(r_,2)
                r(3)=r_(i);
                [u_opt, thetaV, gammaV, X_nfm,  U_nfm_notr, theta] = solve_nfm_over_time(o, X_target, T, thetaV_0, o.cap*ones(1,T), r);
                c_infr_opt=sum(o.c' * U_nfm_notr);
                c_sla_opt= sum(pos( X_sla(:,2:T)-X_nfm(:,2:T)));
                c_trsh_opt= sum(sum(sum(abs(u_opt(:,:,:)))));
                disp(sprintf('\n c_sla_opt=%5.4f, c_infr_opt=%5.4f, c_trsh_opt=%5.4f',...
                    c_sla_opt,  c_infr_opt,  c_trsh_opt));
                % optcost =  r(1)*c_sla_opt + r(2) * c_infr_opt + r(3)*c_trsh_opt
                
                cost_infr = [cost_infr   c_infr_opt];
                cost_trsh =[cost_trsh  c_trsh_opt];
                % fprintf('.');
            end
            [cost_infr' cost_trsh']
            h=figure;
            style  = {'b-.','g--','r-'};
            plot( cost_infr, cost_trsh, style{1}, 'LineWidth',2  );
            xlabel('Infrastructure Cost'); ylabel('Trashing Cost');
            %xlabel('$cost_{infr}$','Interpreter','LaTex')
            %ylabel( '$cost_{trsh}$','Interpreter','LaTex' );
            % [cost_U(2) cost_X(2)]
            for  j=1:3:r_size %   j=[1,round(r_size/2),r_size] %90:3:100
                text( cost_infr(j), cost_trsh(j), sprintf('\\leftarrow r_{trsh}=%4.3f',r_(j)),  'HorizontalAlignment','left')
            end
            title('Tradeoff Curve');
            
            UtilityLib.print_figure(h,9,7,sprintf('figure/tradeoff_curveT%d_rsize%d', T,r_size ));
            
        end
        
        
        function test_mpc_coverage(o)
            T=100;
            %  o.setup_services();
            % N=[250 100 200];
            Z=[1]';
            r=3;
            RT_sla=[0.146]';
            N=200+[1:T]+50*sin([1:5:5*T]/10);
            plot(N);
            % try single norm
            cvx_begin  quiet
            cvx_precision best
            variables  Nopt(T)
            Nopt>=N';
            minimize  r* sum(abs(Nopt(1:T-1) - Nopt(2:T) ))+...  % cost of trashing
                sum(pos(Nopt))%+...     % cos of infrastructure
            % 10*sum(pos(N'-Nopt) )   % sla
            cvx_end
            Nadi=Nopt;
            
            % %  this one ignores the fact that hill is getting deeper and deeper
            nsteps = 30;
            T=30;
            % Xall = zeros(n,nsteps); Uall = zeros(m,nsteps);
            Nall = zeros(1,nsteps);
            for t=1:nsteps
                n=N(t:t+T-1);
                cvx_begin  quiet
                cvx_precision best
                variables  Nopt(T)
                % Nopt(1)==n(1);
                Nopt(1:T)>=n';
                minimize  r* sum(abs(Nopt(1:T-1) - Nopt(2:T) ))+...  % cost of trashing
                    sum(pos(Nopt)) %+...     % cost of infrastructure
                cvx_end
                Nall(t)=Nopt(2);
                Nopt(2)
            end
            %              plot(1:nsteps, N(1:nsteps),'b-', 1:nsteps, Nadi(1:nsteps) ,'r--');
            plot(1:nsteps, N(1:nsteps),'b-', 1:nsteps, Nadi(1:nsteps) ,'r--',1:nsteps, Nall(1:nsteps),'g--');
        end
        
         
        % o= assignment_test_flow_surr()
        % o.test_nfm_mpc_static_compare
        function test_nfm_mpc_static_compare_old(o)
            nsteps = 1;
            T=3;
            o.setup_hosts();
            o.setup_services();
            N=[];
            N=[250 100 200
                100   250 50];
            Z=[1 1]';
            RT_sla=[0.146  0.267]';
            alphaV_0 = zeros(o.num_host , o.num_container);
            
            f_sla = calculate_throuput_objective(o, N, Z, RT_sla, T);
            
            [u, alphaV,betaV,gammaV,~,~] = solve_nfm_over_time(o,f_sla, T, alphaV_0, o.cap*ones(1,T+1));
            
            
            % alphaV
            placement = (alphaV(:,:,:)~=0);
            o.placement = placement;
            % placement
            pl_delta = placement(:,:,2:T+1) -  placement(:,:,1:T)
            o.pl_delta = pl_delta;
            fprintf('initial_deployment=%d \t additions=%d \t removes=%d \t changes=%d \n',...
                sum(sum(sum(pl_delta(:,:,1)))),...
                sum(sum(sum(pl_delta>0))),...
                sum(sum(sum(pl_delta<0))),...
                sum(sum(sum(pl_delta~=0))));
            
            % additions
            [i j k]=ind2sub(size(pl_delta),find(pl_delta(:,:,2:end)>0));
            additions=[i j k+1]
            % removes
            [i j k]=ind2sub(size(pl_delta),find(pl_delta(:,:,2:end)<0));
            removals=[i j k+1]
            o.print_delta_sequece(pl_delta)
            o.print_mat( pl_delta(:,:,1) );
        end
        
         function   test_nfm_mpc_static_compare(o,r)
            nsteps = 1;
            T=4;
            o.setup_hosts();
            o.setup_services();
            N=[];
            N=[0 220 150 200
                 0 100   170 60];
            o.print_mat(N); 
            Z=[1 1]';
            RT_sla=[0.146  0.267]';
            
            %.....
          %  r=[5 1 10];
            % x0=calculate_throuput_objective(o, N(:,1), Z, RT_sla, 1);  
            RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',sum(10)));
            %o.c=rand(o.num_host,o.num_service)+10;
           % o.cost=round(rand(o.num_host,1)*10)+1;
            o.c=rand(o.num_host,o.num_service);
            
           % [theta,thetaV]=o.solve_nfm(x0,r);
           % thetaV_0 =thetaV;  
            thetaV_0 = zeros(o.num_host , o.num_service);
            X_target = calculate_throuput_objective(o, N, Z, RT_sla, T);
            [u, thetaV, gammaV, X_nfm,  U_nfm] = solve_nfm_over_time(o,...
                X_target, T, thetaV_0, o.cap*ones(1,T),r, 'square');
           X_target
            X_nfm
            U_nfm
            bsxfun(@times, o.cap, o.speed_ratio)  
            placement = (thetaV(:,:,1:T)~=0);
            o.placement = placement;
 
            % placement
            pl_delta = placement(:,:,2:T) -  placement(:,:,1:T-1) 
            o.pl_delta = pl_delta;
            for t1=2:T
                 o.draw_service_callgraph_nocontainer(N, Z, RT_sla, true, t1); 
                  waitforbuttonpress
            end 
            
            fprintf('initial_deployment=%d \t additions=%d \t removes=%d \t changes=%d \n',...
                sum(sum(sum(pl_delta(:,:,1)))),...
                sum(sum(sum(pl_delta>0))),...
                sum(sum(sum(pl_delta<0))),...
                sum(sum(sum(pl_delta~=0))));
            
            % additions
            [i j k]=ind2sub(size(pl_delta),find(pl_delta(:,:,2:end)>0));
            additions=[i j k+1]
            % removes
            [i j k]=ind2sub(size(pl_delta),find(pl_delta(:,:,2:end)<0));
            removals=[i j k+1]
            o.print_delta_sequece(pl_delta)
            % o.print_mat( pl_delta(:,:,1) );
        end
        
        function main1(o)
            num_host=3;
            num_service=4;
            num_classes = 3;
            
            RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',sum(10)));
            % rng(sum(100*clock),'v4')
            c=rand(num_host,num_service)*10;
            %c=round(rand(num_host,1)*10);
            %c=round(ones(num_host,1)*10);
            
            % c_r = reshape(c(:,:),[num_host*num_service,1]);
            w=workload();
            nsteps = 1200;
            T=1200;
            tmp= w.get_workload('data/day42_per_min.txt', 1, 23*60)';
            N=[];
            N(1,:) = [tmp(1:nsteps)];
            % lambda(nsteps+1:end)=0;
            N(2,:) =  [sin((1:nsteps)/300)*200+300];
            N(3,:) = [sin(((1:nsteps)+400)/300)*100+400];
            plot(N')
            RT_sla=[10 2 7]';
            Z=[15 15 15]';
            
            % I assume perfect knowlege of workload
            % in computing 'next step desired' throughput
            f_sla = [zeros(3,1) N./((RT_sla + Z ) * ones(1,T))];
            
            % f_sla=[10 2 7]'*ones(1,T) ;
            cap=[60 90 120]';
            
            r1 = [10000 10000 10000];
            
            d= rand(num_service,num_classes);
        end  % main1
        
          function simulate(o)
            % use the actual call graph to build the lqm
          end
    
    end
    
end


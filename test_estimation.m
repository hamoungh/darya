% solve a model over time using stupid linear approach
% try least square, em, mhe, filter and see the accuracy and convergence criteria
classdef test_estimation  < handle 
    
    properties
             num_host
        num_container
        num_service
        num_classes  % 3 4 3 
    % num_inst_types
       
        beta_deployed
        call_graph
        c 
        cost  
        d
        
        N
        Z
        RT_sla  
        cap
        speed_ratio
     
         placement 
         pl_delta
         
         path_str
    end
    
    methods
        function o=test_estimation()     
            o@handle();      
            % http://www.mathworks.com/matlabcentral/answers/12422-macosx-encoding-problem
            cur_file_path=mfilename('fullpath');   [path_str, file_name, ext] = fileparts(cur_file_path);
            o.path_str=path_str; 
             import opera.*;
            import opera.KalmanFilter.*;
        end
        
        function [u, thetaV,gammaV, f_c, U_host_, theta ]=solve_nfm_over_time(o, f_sla, T, thetaV_0, cap, r)
            d=o.d';
            r1 = 1000*ones(1,o.num_classes);
            cvx_begin  quiet  
            cvx_precision best
            variables  thetaV(o.num_host , o.num_service, T+1) ... % state variable
                u(o.num_host , o.num_service, T) ...    % input
                gammaV(o.num_service , o.num_classes,T) ... % gammaV is directly controlable, its input variable
                f_c(o.num_classes,T)...  %  output variable in MPC  based on input
                U_host(o.num_host,T) ... % output, utilization
                consumption_cost(o.num_host,T) ... % cost of output
                consumption_cost_tot(1,T)       % cost of output
            %-------------- constraints ----------------------
            thetaV>=0
            U_host(:,2:T)  <= .8 .* cap(:,2:T) ;  %  o.cap*ones(1,T+1) <= 0 ;
            f_c(:,1:T) >= f_sla;
            %--------------- system dynamics ---------------
            thetaV(:,:,1) ==  thetaV_0
            thetaV(:,:,2:T+1) ==  thetaV(:,:,1:T) + u(:,:,1:T)
            %------------------- system output -----------------
            permute(sum(thetaV(:,:,1:T),1),[2,3,1])  == permute(sum(gammaV, 2),[1,3,2])
            gammaV(:,:,1:T) == ...
                repmat(d,[1,1,(T)]) .* reshape(...
                (ones(o.num_service,1)* reshape(f_c(:,1:T),[1,o.num_classes*(T)])), ...
                [o.num_service,o.num_classes, (T)])
            %--------------------- cost -------------------
            U_host == permute(sum(thetaV(:,:, 1:T),2),[1,3,2]) ;
            consumption_cost_tot ==  o.c'* U_host
            
            minimize(... % time
                r(1) *  sum(pos( f_sla-f_c(:,1:T)))   +  ... % sla violation cost
                r(2) *  sum(consumption_cost_tot(:,1:T)) + ... % infrastructure cost
                r(3) *  sum(sum(sum(sum(abs(u(:,:,:)))))) )% r(3) * sum(sum(sum(abs(u(:,:,:))))) ... % change of control cost summed over host and container
            
            cvx_end
            U_host_=U_host;
            theta= bsxfun(@rdivide, thetaV, sum(thetaV,1));
        end % solve_nfm_over_time
        
           function [f_,rt_,util, d_csh]=solve_lqm(o,theta,N,Z)
            try
                model=OpModel();

                % one delay center for each class 
                model.nodes=[];  model.services=[];  % nodesStr = ['ClientH ']; 
                 scWorkloads=containers.Map(); totN=0; 
                 for c=1:o.num_classes       
                    % nodesStr = [nodesStr,   nodesStr];  
                    model.nodes= [model.nodes, OpNode(sprintf('H_c%d',c),'client',1,1)];  
                    model.services=[model.services, OpService(sprintf('S_c%d',c) , sprintf('T_c%d',c))];  
                    model.containers=[model.containers, OpContainer(sprintf('T_c%d',c) , 1000, sprintf('H_c%d',c) , 'false' )];     
                    % scWorkloads=containers.Map({sprintf('class%d',c)},{OpClosedWorkload(N,Z)});
                    roundN=round(N(c));
                    scWorkloads(sprintf('class%d',c)) = OpClosedWorkload(roundN,Z(c)); 
                    totN = totN + roundN; 
                    model.scenarios= [model.scenarios, OpScenario(sprintf('class%d',c) , sprintf('S_c%d',c))];  
                 end   
                model.workload = OpWorkload(totN,scWorkloads);  
                    
                for i=1:o.num_host
                    model.nodes=[model.nodes,  OpNode(sprintf('H%d',i),'server', o.cap(i) ,1)];
                end
                
                % for now i put each service on a container                        
                d_csh= permute(  repmat(o.d, [ 1,1 , o.num_host] ) , [1,2,3])  .* ... 
                                 permute( repmat(theta, [1 , 1, o.num_classes]), [3,2,1]  ); 
            
                for s=1:o.num_service               
                     for h=1:o.num_host 
                        if (theta(h,s)~=0)       
                            model.services=[model.services,  OpService(sprintf('S%d_%d',s,h), sprintf('T%d_%d',s,h))];   
                            model.containers=[model.containers,   OpContainer(sprintf('T%d_%d',s,h), 1000, sprintf('H%d',h), 'true' )];
                            for c=1:o.num_classes      
                                         scs=model.scenarios; 
                                         sc=scs(c);
                                         sc.addCall( OpCall(sprintf('S_c%d',c) , sprintf('S%d_%d',s,h), 1, d_csh(c,s,h), 0)  ); 
                             end  
                         end 
                    end
                end

             model.clusters=[OpCluster('ClientCluster',[model.containers(1:end)])]; 
                 
                nodesStr = [];
                for i=1:o.num_host+o.num_classes   
                    nodesStr = [nodesStr, model.nodes(i).name, ' ']; 
                end
                model.networks = OpNetwork(nodesStr);
            catch err
                rethrow(err);
            end

            %+++++++++++++++++
             model.solve();
            util = [];
            for i=o.num_classes+1:length(model.nodes)
                util=[util model.nodes(i).cpuUtilization];
            end
      ;
            
            f_=[]; rt_=[]; 
            for c=1:o.num_classes       
                f_=[f_ ; model.scenarios(c).throughput]; 
                rt_=[rt_; model.scenarios(c).responseTime];   
            end 
            
        end % solve_lqm
       
         function [ X_nfm, RT_nfm, U] = calculate_nfm(o, thetaV, gammaV, D, N, Z)  % linearize output to state relation knowing the current state thetaV, gammaV 
            f_sla = [N./((RT_sla + Z ) * ones(1,T))];   
            n=x(z+r)
             d=o.d';    
             r1 = 1000*ones(1,o.num_classes); 
            cvx_begin quiet
            variables gammaV(o.num_service,o.num_classes) ...
                X_c(o.num_classes,1) 
                %thetaV(o.num_host,o.num_service)...
   
            sum(thetaV,1)' == sum(gammaV ,2) 
             gammaV == d .* (ones(o.num_service,1)*X_c')  
             
             U_host == permute(sum(thetaV(:,:, 1:T),2),[1,3,2]) ;
             %   X_c >= X_target       
            minimize sum( sum(o.c .* thetaV,2)) 
            cvx_end  
            theta= bsxfun(@rdivide, thetaMu, sum(thetaMu,1));
        end % solve_nfm
       
       
        function calculate_jacobian(o) 
            
        end 
        
        function  Y=calculate_contact_level(o,call) % write using optimization
       
            con=zeros(o.num_service+o.num_classes,o.num_service+o.num_classes);
            for a=1:o.num_service
                con_old=con;
                con= con+ call^a;
                if (con_old==con)
                    break;
                end;
            end
            Y=con(1:o.num_classes, o.num_classes+1:o.num_classes+o.num_service);
        end
        
           function f_sla = calculate_throuput_objective(o, N, Z, RT_sla, T) 
                 % I assume perfect knowlege of workload 
                 % in computing 'next step desired' throughput 
                  f_sla = [N./((RT_sla + Z ) * ones(1,T))];   
           end
           
           function RT = convert_X_to_RT(o, N, Z, X, T)
                 % I assume perfect knowlege of workload 
                 % in computing 'next step desired' throughput 
                   RT = N./X - Z* ones(1,T); 
         end
           
        function [optcost_tr, mpccost, optcost] = test_mpc_overal(o, nsteps , r_trsh_cost)
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
            T=nsteps; 
  
            thetaV_0 =  o.cap;  %[0   8.0000   0.0360  0]';
           % N=200+50*sin([1:1*(nsteps)]/50)+...
          %  N =   cumsum( 5*randn(o.num_classes,T)')';     % +[1:nsteps] 
          %  N = 200+N-min(N);   
          %  X_target = calculate_throuput_objective(o, N, Z, RT_sla, nsteps);
          X_target =   cumsum( 5*randn(o.num_classes,T)')';  
          X_target = 200+X_target-min(X_target);     
          r=[50 1 r_trsh_cost]; 
         [u_opt, thetaV, gammaV, X_nfm,  U_nfm_notr, theta] = solve_nfm_over_time(o, X_target, T, thetaV_0, o.cap*ones(1,T), r);
         Z=[1]';
         RT =[0.146]';
         N = X_target .* ((Z+RT)* ones(1,T));
            
%       %     H=calculate_jacobian(o);   % change of throughput, response time, utilization with respect to D, N, Z 
%            
%            subplot(4,1,1);  plot(N); title('N','Interpreter','none'); 
%            subplot(4,1,2);  plot(U_nfm_notr'); title('U_nfm_notr','Interpreter','none');  % util of each server 
%            subplot(4,1,3);  plot(X_nfm(2:end)); title('X_nfm','Interpreter','none');   % throughput  
%            RT = convert_X_to_RT(o, N, Z, X_nfm, T); % response time 
%            subplot(4,1,4);   plot(RT(2:end)); title('RT')
%             
           % http://www.mathworks.com/matlabcentral/answers/12422-macosx-encoding-problem
           import opera.*;
           import opera.KalmanFilter.*;
           theModel =  OperaModel(); % OperaModel
           theModel.setModel(strcat(o.path_str,'/opera/output/input.pxl'));
           kalmanConfig = KalmanConfiguration(); %KalmanConfiguration
           kalmanConfig.withConfigFile(strcat(o.path_str,'/opera/output/mykalman.config'))...
                .withSetting(KalmdanConfiguration.ITERATIONS_MAX, '5')...
                .withModel(theModel)...
                .withSetting(KalmanConfiguration.FILE_TRACE, strcat(o.path_str,'/opera/output/mykalman.trace.txt'));
                %     %    .withSetting(KalmanConfiguration.FILE_MODEL_RESULTS, strcat(o.path_str,'/opera/output/output.xml'))...
           theEstimator =  KalmanEstimator(kalmanConfig); %KalmanEstimator
                 
           
%      for loop
%          % get new measurement 
%           [X_,rt_,util]=solve_lqm(o,theta(:,:,2),N(2),Z);  
%           % put the mesurements in the kalman as defaults in the config file  
%           %
%           % estimate through kalman, it will update its internal model u  passed to 
%          estimationResults = theEstimator.EstimateModelParameters();
%          % get the estimated outputs 
%          estimationResults .metricsMeasured 
%          % get the estimated parameters 
%          % 
%          % calc overall error , if its too high we need more clusters 
%          
%          % calc the estimation error, if its too high and doesnt converge we need less clusters 
%      end 
%          
          % update measurements 
           
%             c_infr_opt=sum(o.c' * U_nfm_notr);
%             c_sla_opt= sum(pos( X_sla(:,2:T)-X_nfm(:,2:T)));
%             c_trsh_opt= sum(sum(sum(abs(u_opt(:,:,:)))));
%             disp(sprintf('\n c_sla_opt=%5.4f, c_infr_opt=%5.4f, c_trsh_opt=%5.4f',...
%                 c_sla_opt,  c_infr_opt,  c_trsh_opt));
%             optcost =  r(1)*c_sla_opt + r(2) * c_infr_opt + r(3)*c_trsh_opt
        end % test
        
        function test_simple(o,nsteps)
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
            T=nsteps; 
  
            thetaV_0 =  o.cap;  %[0   8.0000   0.0360  0]';
           % N=200+50*sin([1:1*(nsteps)]/50)+...
          %  N =   cumsum( 5*randn(o.num_classes,T)')';     % +[1:nsteps] 
          %  N = 200+N-min(N);   
          %  X_target = calculate_throuput_objective(o, N, Z, RT_sla, nsteps);
          X_target =   cumsum( 5*randn(o.num_classes,T)')';  
          X_target = 200+X_target-min(X_target);     
    %      [u_opt, thetaV, gammaV, X_nfm,  U_nfm_notr, theta] = solve_nfm_over_time(o, X_target, T, thetaV_0, o.cap*ones(1,T), r);
         Z=[1]';
         RT =[0.146]';
         N = X_target .* ((Z+RT)* ones(1,T));
         theta = [0.0000 0.5898 0.4102   0.0000]'; 
    
          %  for loop
         % get new measurement 
          [X_,rt_,util, d_csh]=solve_lqm(o,theta,N(2),Z);  
          % put the mesurements in the kalman as defaults in the config file  
          path=KalmanConfigRenderer().render_and_write(X_,rt_,util, d_csh);
          % estimate through kalman, it will update its internal model u  passed to 
      %   estimationResults = theEstimator.EstimateModelParameters();
         % get the estimated outputs 
      %   estimationResults .metricsMeasured 
         % get the estimated parameters 
         % 
         % calc overall error , if its too high we need more clusters 
         
         % calc the estimation error, if its too high and doesnt converge we need less clusters 
        %  end  % loop 
  
        end %test......
        
    end  %methods
    
    end % class 
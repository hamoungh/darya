
classdef autoscaling_mpc < handle
    properties  
        num_clusters
        num_service
        num_classes  
        num_host
        beta_deployed 
        cap  
    end 
    
    methods
         function f_sla = calculate_throuput_objective(o, N, Z, RT_sla, T) 
                 % I assume perfect knowlege of workload 
                 % in computing 'next step desired' throughput 
                  f_sla = [N./((RT_sla + Z ) * ones(1,T))];   
         end
         
           function [f_,rt_]=solve_lqm(o,d,theta,N,Z)
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
                d_csh= permute(  repmat(d, [ 1,1 , o.num_host] ) , [1,2,3])  .* ... 
                                 permute( repmat(theta, [1 , 1, o.num_classes]), [3,2,1]  ); 
                 % disp('demand'); 
                 % didsp(permute(d_csh,[1,3,2])  ); 
                % this is going to be more complex for different examples,
                % later on, maybe comes from a matrix
                % but now its just equal calls from the client to all services
                for s=1:o.num_service               
                     for h=1:o.num_host 
                        if (theta(h,s)~=0)       
                            model.services=[model.services,  OpService(sprintf('S%d_%d',s,h), sprintf('T%d_%d',s,h))];   
                            model.containers=[model.containers,   OpContainer(sprintf('T%d_%d',s,h), 1000, sprintf('H%d',h), 'true' )];
                            for c=1:o.num_classes      
                                    % OpCall(caller, callee, invocations, CPUDemand, DiskDemand)       
                                         scs=model.scenarios; 
                                         sc=scs(c);
                                         sc.addCall( OpCall(sprintf('S_c%d',c) , sprintf('S%d_%d',s,h), 1, d_csh(c,s,h), 0)  ); 
                                         %model.calls=[model.calls...
                                          %OpCall('ClientS', sprintf('S%d_%d',s,h), 1*theta(h,s), o.d(s,c), 0)];     
                                         % ];           
                             end  
                         end 
                    end
                end
                
                %  model.clusters=[OpCluster('ClientCluster',[model.containers(1)]),...
                %  OpCluster('TaskCluster',[model.containers(2:end)])];
             model.clusters=[OpCluster('ClientCluster',[model.containers(1:end)])]; 
                 
                nodesStr = [];
                for i=1:o.num_host+o.num_classes   
                    % nodesStr = [nodesStr, strcat(model.nodes(i).name,{' '})]; 
                    nodesStr = [nodesStr, model.nodes(i).name, ' ']; 
                    % nodes=[nodes sprintf('H%d ',i)];
                end
                model.networks = OpNetwork(nodesStr);
            catch err
                rethrow(err);
            end
       
            
            %+++++++++++++++++
            model.solve();
            %         util = [];
            %         for i=3:length(opM.nodes)
            %             util=[util opM.nodes(i).cpuUtilization];
            %         end
            %         util_ = mean(util);
            f_=[]; rt_=[]; 
            for c=1:o.num_classes       
                f_=[f_ ; model.scenarios(c).throughput]; 
                rt_=[rt_; model.scenarios(c).responseTime];   
            end 
            
        end % solve_lqm
      
          function [f_,rt_]=solve_lqm3(o,d, betaV,N,Z) 
            multi=ceil(sum(betaV,2));
            % multi=sum(betaV,2); 
            beta = bsxfun(@rdivide, betaV, sum(betaV,1));  
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
                    
                for cl=1:o.num_clusters 
                    model.nodes=[model.nodes,  OpNode(sprintf('H%d',cl),'server', multi(cl) ,1)];
                end
                
                % for now i put each service on a container                        
                d_csh= permute(  repmat(d, [ 1,1 , o.num_clusters] ) , [1,2,3])  .* ... 
                                 permute( repmat(beta, [1 , 1, o.num_classes]), [3,2,1]  ); 
                 % disp('demand'); 
                 % didsp(permute(d_csh,[1,3,2])  ); 
                % this is going to be more complex for different examples,
                % later on, maybe comes from a matrix
                % but now its just equal calls from the client to all services
                for s=1:o.num_service               
                     for h=1:o.num_clusters
                        if (beta(h,s)~=0)       
                            model.services=[model.services,  OpService(sprintf('S%d_%d',s,h), sprintf('T%d_%d',s,h))];   
                            model.containers=[model.containers,   OpContainer(sprintf('T%d_%d',s,h), 1000, sprintf('H%d',h), 'true' )];
                            for c=1:o.num_classes      
                                    % OpCall(caller, callee, invocations, CPUDemand, DiskDemand)       
                                         scs=model.scenarios; 
                                         sc=scs(c);
                                         sc.addCall( OpCall(sprintf('S_c%d',c) , sprintf('S%d_%d',s,h), 1, d_csh(c,s,h), 0)  ); 
                                         %model.calls=[model.calls...
                                          %OpCall('ClientS', sprintf('S%d_%d',s,h), 1*theta(h,s), o.d(s,c), 0)];     
                                         % ];           
                             end  
                         end 
                    end
                end
                
                %  model.clusters=[OpCluster('ClientCluster',[model.containers(1)]),...
                %  OpCluster('TaskCluster',[model.containers(2:end)])];
             model.clusters=[OpCluster('ClientCluster',[model.containers(1:end)])]; 
                 
                nodesStr = [];
                for n=1:o.num_clusters+o.num_classes   
                    % nodesStr = [nodesStr, strcat(model.nodes(i).name,{' '})]; 
                    nodesStr = [nodesStr, model.nodes(n).name, ' ']; 
                    % nodes=[nodes sprintf('H%d ',i)];
                end
                model.networks = OpNetwork(nodesStr);
            catch err
                rethrow(err);
            end
       
            
            %+++++++++++++++++
            model.solve();
            %         util = [];
            %         for i=3:length(opM.nodes)
            %             util=[util opM.nodes(i).cpuUtilization];
            %         end
            %         util_ = mean(util);
            f_=[]; rt_=[]; 
            for c=1:o.num_classes       
                f_=[f_ ; model.scenarios(c).throughput]; 
                rt_=[rt_; model.scenarios(c).responseTime];   
            end 
            
        end % solve_lqm
      
        % in that host we are going to host all the
        % services of the cluster with the demand  splited 
        % for each service, all classes demands are splitted
         % Graph Algorithms in the Language of Linear Algebra  
         % Bipartite Graphs and Their Applications, page 15 is about matrix
         % graph mapping 
        % seys what portion of cluster multiplicities are given to services
        % betaV(o.num_clusters,o.num_service)    
        % gammaV(o.num_service,o.num_classes)   
          function [theta,thetaV,cap]=calculate_host_from_cluster (o, machineSize, betaV)
                % compute a theta(o.num_host,o.num_service) based on the betaV(o.num_clusters,o.num_service)    
                host_per_cluster=ceil(sum(betaV,2)./machineSize);
                thetaV=[]; cap=[]; 
                for cl=1:o.num_clusters   
                    for host=1:host_per_cluster(cl) 
                        thetaV=[thetaV; betaV(cl,:)/host_per_cluster(cl) ];      
                        cap=[cap machineSize(cl)]; 
                    end
                end
              theta= bsxfun(@rdivide, thetaV, sum(thetaV,1));
          end  %function 
     
       function [betaV]=solve_nfm(o, d, X_target)  
            % d=o.d';    
           % multiplicity of cluster mult its server multi should support the util on the service 
           % here i assume i know the D_c,s=sum D_c,s,h from estimation 
           % so here we talk about mult of software resource 
           % assume there is no hardware contention and each software
           % replice is placed  on a separate hardware 
            cvx_begin quiet   
                variables  betaV(o.num_clusters,o.num_service)... 
                    gammaV(o.num_service,o.num_classes) ...
                    X_c(o.num_classes,1) 
                X_c >= X_target; 
                betaV>=0;  
                betaV .* o.beta_deployed == betaV;    
                sum(betaV,1)' == sum(gammaV,2)    
                gammaV == d' .* (ones(o.num_service,1)*X_c')  
                 %   X_c >= X_target       
                minimize sum( sum(betaV,2))
            cvx_end  
            % note that here each service is only serviced by one cluster 
            % o.beta_deployed = bsxfun(@rdivide, thetaMu, sum(thetaMu,1))
        end % solve_nfm
  
          function [betaV,X_c,u,inf_cost,tr_cost]=solve_nfm_over_time(o, d, X_target,betaV_0,T,r)  
            cvx_begin    
                variables  betaV(o.num_clusters,o.num_service, T+1)...  % % state variable  
                    gammaV(o.num_service,o.num_classes, T) ...  % % gammaV is directly controlable, its input variable  
                    X_c(o.num_classes,T) ...    %    output variable  based on input  
                    u(o.num_clusters,o.num_service, T)   % input    
                betaV(:,:,1) ==  betaV_0 
                betaV(:,:,2:T+1) ==  betaV(:,:,1:T) + u(:,:,1:T)   
                X_c (:,1:T)>= X_target;  
                betaV>=0;  
                betaV .*  repmat(o.beta_deployed, [1 1 T+1]) == betaV;     
                permute(sum(betaV(:,:,1:T),1),[2,3,1])  == permute(sum(gammaV, 2),[1,3,2])     
                % gammaV == d' .* (ones(o.num_service,1)*X_c')   
                gammaV(:,:,1:T) == ...   
                     repmat(d',[1,1,T]) .* reshape(...
                            (ones(o.num_service,1)* reshape(X_c(:,1:T),[1,o.num_classes*T])), ...
                            [o.num_service,o.num_classes, T])
                        
                minimize r(2)*sum(sum( sum(betaV,2),1),3)...    %    + r(3)*sum(sum(sum(square(u),2),1),3) 
                                 + r(3)*sum(sum(sum(abs(u),2),1),3)   % sum(abs(u),2) are the clusters multiplicity which we call m in the paper 

            cvx_end  
            inf_cost = permute(sum(sum(betaV(:,:, 1:T),2),1), [1,3,2]) ; 
           % tr_cost=permute(sum(sum(square(u),1),2), [1 3 2]); 
            tr_cost=permute(sum(sum(abs(u),1),2), [1 3 2]); 
            % note that here each service is only serviced by one cluster 
            % o.beta_deployed = bsxfun(@rdivide, thetaMu, sum(thetaMu,1))
        end % solve_nfm
        
    function  Y=calculate_contact_level(o,call) % write using optimization  
              % total direct and indirect requests for a service
                % described in 'requests per user response'
                % Yuser1DB2Serv = 1 x 0.3 x (0.7 + 0.3 x 1.4) = 0.336
                % found by following all paths from the user class to the service.
                % calculating the contact level between class and service
                % nodes
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
        
    function test_simple_case(o) 
                o.num_service= 2;
                o.num_clusters= 2; 
                o.num_classes = 1;  
                
                beta=zeros(o.num_clusters,o.num_service);
                beta(1,1)=1; 
                beta(2,1)=1; beta(2,2)=1; 
                o.beta_deployed=beta; 

                call=zeros(o.num_service+o.num_classes,o.num_service+o.num_classes);
                  call=[0 1 1
                           0 0 0
                           0 0 0];  
                Y=calculate_contact_level(o,call);
                
                D =  [0.04 0.02];  
                d= Y .* (ones(o.num_classes,1) * D); 
   
                  N=[250]; 
                  Z=[1]';
                  RT_sla=[0.146]'; 
                  X_target = [N./(RT_sla + Z )];   
              % X_target = calculate_throuput_objective(o, N, Z, RT_sla, T); 
                 [betaV]=solve_nfm(o, d, X_target) 
               % model 1  
                 machineSize =   [9 9]';% ceil(sum(betaV,2));     % 
                 machineSize= ones(o.num_clusters,1);  
                  [theta,thetaV,cap]=calculate_host_from_cluster (o, machineSize, betaV); 
                  o.cap=cap; 
                  o.num_host=size(theta,1); 
                  [f_,rt_]=solve_lqm(o, d, theta,N,Z)   
                 % model 2  
                  [f_,rt_]=solve_lqm3(o, d, machineSize, betaV,N,Z)   
                 
    end %test_simple_case(o) 
    
     function [X_lqm,r_lqm]=eveluate_lqm(o, d, betaV , N,  Z)  
             X_lqm=[]; r_lqm=[]; 
              for t=1:size(N,2)  
                     [x_q, r_q] = solve_lqm3(o, d, betaV(:,:,t),N(:,t) ,Z); 
                    X_lqm=[x_q X_lqm]; 
                    r_lqm=[r_q r_lqm]; 
                   %  disp(sprintf('solving lqm for t=%d',t));  
                    fprintf('.');    
              end
     end
         
      function test_simple_case_over_time(o,nsteps, r_trsh_cost)     
                o.num_service= 2;
                o.num_clusters= 2; 
                o.num_classes = 1;  
                
                beta=zeros(o.num_clusters,o.num_service);
                beta(1,1)=1; 
                beta(2,1)=1; beta(2,2)=1; 
                o.beta_deployed=beta; 

                call=zeros(o.num_service+o.num_classes,o.num_service+o.num_classes);
                call=[0 1 1
                         0 0 0
                         0 0 0];  
                Y=calculate_contact_level(o,call);
                
                D =  [0.04 0.02];  
                d= Y .* (ones(o.num_classes,1) * D);   
                  Z=[1]';
                  RT_sla=[0.146]'; 
                % nsteps=200; 
                 T=nsteps; 
                  betaV_0 =  [4.3630         0
                                      4.3630    4.3630];  % [0   8.0000   0.0360  0]';    

                     tmp= workload().get_workload('data/day42_per_min.txt', 1, 23*60)'; 
                     % tmp = workload().get_workload('data/day66_per_min.txt', 1, 23*60);
                     N = [tmp(1:nsteps)];
                     a = 1;
                     b = [1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8];
                     N_smooth = filter(b,a,N);
                    % plot(N_smooth); hold on; 
                 %    N=200+[1:nsteps]+50*sin([1:5:5*(nsteps)]/10); plot(N); 
                  X_target = calculate_throuput_objective(o, N_smooth, Z, RT_sla, nsteps); 
              
                  % the optimal case 
                  r=[50 1 r_trsh_cost];    
                  [betaV_opt,X_nfm_opt,u_opt, inf_cost_opt,tr_cost_opt] = solve_nfm_over_time(o, d, X_target , betaV_0,nsteps, r); 
                  plot(ceil(permute(sum(betaV_opt,2),[3,1,2])));
                   title('#Servers'); legend('tomcat cluster','tomcat+apache cluster' ); 
                   hold on; 
                  figure; plot(N); 
                  
                c_infr_opt=sum(sum(sum(betaV_opt,2),1),3); 
                c_sla_opt= sum(pos( X_target(:,2:T)-X_nfm_opt(:,2:T))); 
                c_trsh_opt= sum(sum(sum(abs(u_opt(:,:,:)))));
                disp(sprintf('\n c_sla_opt=%5.4f, c_infr_opt=%5.4f, c_trsh_opt=%5.4f',...
                                c_sla_opt,  c_infr_opt,  c_trsh_opt));  
               optcost =  r(1)*c_sla_opt + r(2) * c_infr_opt + r(3)*c_trsh_opt
               
%                % the case that the trashing is ignored 
%                     r_trsh=[50 1 0.01];
%                   [betaV_tr,X_nfm_tr,u_opt_tr, inf_cost_tr,tr_cost_tr] = solve_nfm_over_time(o, d, X_target , betaV_0,nsteps, r_trsh); 
%               %   [X_lqm,r_lqm]=eveluate_lqm(o, d, betaV , N,  Z) 
%                
%              % compare cost of infrastructure 
%             h=figure; 
%             plot(1:T, inf_cost_tr  ,'-bs',...
%                     1:T, inf_cost_opt  ,'-ro');  
%             hleg1 = legend('trashing ignored','optimal' ); 
%             title('Cost of infrastructure'); 
%             xlabel('Time step'); ylabel('Infrastructure Cost');
%             UtilityLib.print_figure(h,9,7,sprintf('figure/ae_infrastructure_cost_nsteps%d_T%d_rTRSH%d', nsteps, nsteps,r_trsh_cost));
%            
%             % compare cost of trashing 
%               h=figure; 
%               plot(1:T-1, tr_cost_tr(1:T-1),  '-bs',...
%                   1:T-1,  tr_cost_opt(1:T-1)  ,'-ro'); 
%                hleg1 = legend('trashing ignored','optimal' ,'mpc');  
%                title('Cost of trashing'); 
%                xlabel('Time step'); ylabel('Trashing Cost');
%               UtilityLib.print_figure(h,9,7,sprintf('figure/ae_trashing_cost_nsteps%d_T%d_rTRSH%d', nsteps, nsteps,r_trsh_cost));  
    end %test_simple_case(o) 
  
     function test_simple_case_over_time_clusters(o,nsteps )     
                o.num_service= 2;
                o.num_clusters= 2; 
                o.num_classes = 1;  
                
                beta=zeros(o.num_clusters,o.num_service);
                beta(1,1)=1; 
                beta(2,1)=1; beta(2,2)=1; 
                o.beta_deployed=beta; 

                call=zeros(o.num_service+o.num_classes,o.num_service+o.num_classes);
                call=[0 1 1
                         0 0 0
                         0 0 0];  
                Y=calculate_contact_level(o,call);
                
                D =  [0.04 0.02];  
                d= Y .* (ones(o.num_classes,1) * D);   
                  Z=[1]';
                  RT_sla=[0.146]'; 
                % nsteps=200; 
                 T=nsteps; 
                  betaV_0 =  [4.3630         0
                                      4.3630    4.3630];  % [0   8.0000   0.0360  0]';    

                     tmp= workload().get_workload('data/day42_per_min.txt', 1, 23*60)'; 
                     % tmp = workload().get_workload('data/day66_per_min.txt', 1, 23*60);
                     N = [tmp(1:nsteps)];
                     a = 1;
                     b = [1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8];
                     N_smooth = filter(b,a,N);
                    % plot(N_smooth); hold on; 
                 %    N=200+[1:nsteps]+50*sin([1:5:5*(nsteps)]/10); plot(N); 
                  X_target = calculate_throuput_objective(o, N_smooth, Z, RT_sla, nsteps); 
              
                  % the optimal case 
                  r=[50 1 4];    
                  [betaV_opt,X_nfm_opt,u_opt, inf_cost_opt,tr_cost_opt] = solve_nfm_over_time(o, d, X_target , betaV_0,nsteps, r); 
                  subplot(3,1,1); plot(ceil(permute(sum(betaV_opt,2),[3,1,2]))); 
                  title('#Servers, r_2=4'); legend('tomcat cluster','tomcat+apache cluster' ); 
                  
                  r=[50 1 40];    
                  [betaV_opt,X_nfm_opt,u_opt, inf_cost_opt,tr_cost_opt] = solve_nfm_over_time(o, d, X_target , betaV_0,nsteps, r); 
                  subplot(3,1,2); plot(ceil(permute(sum(betaV_opt,2),[3,1,2]))); 
                   title('#Servers, r_2=40'); legend('tomcat cluster','tomcat+apache cluster' ); 
                   
                   subplot(3,1,3);  plot(N); 
                    title('Number of Users'); 
     end 
     
       function draw_tradeoff_curve(o) 
                  o.num_service= 2;
                o.num_clusters= 2; 
                o.num_classes = 1;  
                
                beta=zeros(o.num_clusters,o.num_service);
                beta(1,1)=1; 
                beta(2,1)=1; beta(2,2)=1; 
                o.beta_deployed=beta; 

                call=zeros(o.num_service+o.num_classes,o.num_service+o.num_classes);
                call=[0 1 1
                         0 0 0
                         0 0 0];  
                Y=calculate_contact_level(o,call);
                
                D =  [0.04 0.02];  
                d= Y .* (ones(o.num_classes,1) * D); 
                  Z=[1]';
                  RT_sla=[0.146]'; 
                 nsteps=5; 
                  betaV_0 =  [4.3630         0
                                      4.3630    4.3630];  % [0   8.0000   0.0360  0]';    
                  N=200+[1:nsteps]+50*sin([1:5:5*(nsteps)]/10); 
                  X_sla = calculate_throuput_objective(o, N, Z, RT_sla, nsteps);                  
                  X_target=X_sla;  

                cost_infr=[]; cost_trsh=[];  
                r_size=40;
                % r_=logspace( 0, 2, r_size ) 
                 r_=linspace(1,100,r_size); 
                r=[50 1 0]; 
             
               
             
                
                for i=1:size(r_,2)
                    r(3)=r_(i); 
                   %     [betaV]=solve_nfm_over_time(o, d, X_target , betaV_0,nsteps)
                    % [u_opt, thetaV, gammaV, X_nfm,  U_nfm_notr, theta] = solve_nfm_over_time(o, X_target, T, thetaV_0, o.cap*ones(1,T), r);  
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

         end % function 
    end %methods 
end 
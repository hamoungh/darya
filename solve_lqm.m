function [f_,rt_]=solve_lqm(o,theta,N,Z)
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
classdef testmpc1 < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
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
        c_host_sparse
        c_all_sparse
        d
        
        N
        Z
        RT_sla
        
        %f_sla=[10]';%[10 2 7]';
        % cap becomes the multiplicity of the host later on
        cap
        speed_ratio
        %  r1
        
        % outputs after solving the model
        placement
        pl_delta
    end
    
    methods(Static)
        function test_java()
            disp('hi'); 
        end 
    end
    
    methods
        
        function f_sla = calculate_throuput_objective(o, N, Z, RT_sla, T)
            % I assume perfect knowlege of workload
            % in computing 'next step desired' throughput
            f_sla = [N./((RT_sla + Z ) * ones(1,T))];
        end
        
        function [u, thetaV_round,gammaV, f_c, U_host_, theta, thetaV]=solve_nfm_over_time(...
                o, f_sla, T, thetaV_0, cap, r, objective )
            d=o.d';
            headroom=0.65; 
            cvx_begin       quiet
            cvx_precision high %  low medium default high best
            cvx_solver    sdpt3 %sedumi %
            
            variables  thetaV(o.num_host , o.num_service, T+1) ... % state variable
                u(o.num_host , o.num_service, T) ...    % input
                gammaV(o.num_service , o.num_classes,T) ... % feed forward input variable, directly controlable
                f_c(o.num_classes,T)...  %  output variable in MPC  based on input
                U_host(o.num_host,T) ... % output, utilization
                consumption_cost(o.num_host,T) ... % cost of output
                consumption_cost_tot(1,T)       % cost of output
            
            
            %-------------- constraints ----------------------
            thetaV>=0 %r
            % sum(alphaV,2)<=o.cap
            % permute(sum(thetaV(:,:, 1:T+1),2),[1,3,2])
            % U_host -  cap    <= 0 ;  %  o.cap*ones(1,T+1) <= 0 ;
            U_host -   headroom * bsxfun(@times, cap, o.speed_ratio)   <= 0 ;
            % f_c(:,1:T) >= f_sla;    % have this
            %--------------- system dynamics ---------------
            thetaV(:,:,1) ==  thetaV_0%just added
            thetaV(:,:,2:T+1) ==  thetaV(:,:,1:T) + u(:,:,1:T)
            
            %------------------- system output -----------------
            % betaV .* repmat(o.beta_deployed,[1,1,T+1]) == betaV;
            % gamma state var can be controlled
            permute(sum(thetaV(:,:,1:T),1),[2,3,1])  == permute(sum(gammaV, 2),[1,3,2])
            
            % here im talking about each class node
            % perfect information assumption assumes we know input or workload (N) ahead
            % and proper output (f_sla) can be
            % describing how gamma is going to affect next step's f_c
            % pretty much like 'X(:,2:T+1) == A*X(:,1:T)+B*U' but we are stateless
            %
            % gammaV == d .* (ones(o.num_service,1)*f_c')
            gammaV(:,:,1:T) == ...
                repmat(d,[1,1,(T)]) .* reshape(...
                (ones(o.num_service,1)* reshape(f_c(:,1:T),[1,o.num_classes*(T)])), ...
                [o.num_service,o.num_classes, (T)])
            %   f_c(:,1) == zeros(o.num_classes,1);
            
            %--------------------- cost --------------------
            % minimize sum(sum(o.c.*alphaV)) +  o.r1 * pos( f_sla - f_net)
            % for linear case use pos instead of square_pos
            % not that cost is affine to the resources and flows (for now) so I
            % could extract it but the sla penalty is not (because of pos function) , so I have to put it in
            % objective
            U_host == permute(sum(thetaV(:,:, 1:T),2),[1,3,2]) ;
            %consumption_cost == o.c* U_host %permute(sum( repmat(o.c,[1 1 T]) .*thetaV,2),[1,3,2])
            consumption_cost_tot ==   sum( permute(sum( repmat(o.c(:,:,1), [1,1,T]) .* thetaV(:,:, 1:T),2),[1,3,2]) ,1);
            % consumption_cost_tot ==  sum(consumption_cost,1);
            % minimize( sum(consumption_cost_tot ,2))
            
            %----------------- objective ------------------  make that
            %500 a 0 if you want no cost
            if (strcmpi(objective,'abs'))
                minimize(...
                    r(1) *  sum(sum(pos( f_sla-f_c(:,1:T)),1),2)   +  ... % sla violation cost
                    r(2) *  sum(consumption_cost_tot(:,1:T)) + ... % infrastructure cost
                    r(3) *  sum(sum(sum(abs(u),2),1),3) )
            elseif  (strcmpi(objective,'square'))
                minimize(...
                    r(1) *  sum(sum(pos( f_sla-f_c(:,1:T)),1),2)   +  ... % sla violation cost
                    r(2) *  sum(consumption_cost_tot(:,1:T)) + ... % infrastructure cost
                    r(3) *  sum(sum_square(sum(u,2),1),3) ) % u=o.num_host , o.num_service, T
            elseif  (strcmpi(objective,'no_cost'))
                minimize(...
                    r(1) *  sum(sum(pos( f_sla-f_c(:,1:T)),1),2)   +  ... % sla violation cost
                    r(2) *  sum(consumption_cost_tot(:,1:T)) )
            end
            cvx_end
            
            %               r(1) *  sum(pos( f_sla-f_c(:,1:T)))
            %               r(2) *  sum(consumption_cost_tot(:,1:T))
            %               r(3) *  sum(sum(sum(sum(abs(u(:,:,:))))))
            %   u= round((10^4).*u)/(10^4);
            permute(sum(thetaV(:,:, 1:T),2),[1,3,2]); 
            % thetaV_round used to be thetaV, and it made it really sparse 
            % ( zero values), but now it's just some epsilon above zero
            % thetaV_orig=thetaV;
            thetaV_round = round((10^2).*thetaV)/(10^2); 
            U_host_=permute(sum(thetaV_round,2),[1 3 2]); %U_host_=U_host;
            theta= bsxfun(@rdivide, thetaV_round, sum(thetaV_round,1));
            %  consumption_cost % utilization of hosts
            % mu_host = permute(sum(thetaV(:,:, 1:T+1),2),[1,3,2]) ;
            
                      c_infr=consumption_cost_tot(:,1:T);
                      c_sla= pos( f_sla-f_c(:,1:T));
                      c_trsh= sum(sum(sum(abs(u(:,:,:)))));
            
            
        end % solve_nfm_over_time
        
        function [theta,thetaV,gammaV]=solve_nfm(o, X_target,r, cost)
            d=o.d';
            %e     r1 = 1000*ones(1,o.num_classes);
            cvx_begin quiet
            % f_serv is the requested throughput of the user class
            % d is flow ratio parameters which convert the class flow f_c, in
            % units of user requests/sec, to demand flows ?sc for services
            
            % gamma is service rate in terms of cpu cycles per sec
            %  outputed from services to user classes
            variables  thetaV(o.num_host,o.num_service)...
                gammaV(o.num_service,o.num_classes) ...
                X_c(o.num_classes,1)
            
            thetaV>=0
            sum(thetaV,2)<=o.cap .* o.speed_ratio
            
            sum(thetaV,1)' == sum(gammaV ,2)
            % For each single user request by class c, a demand of d_sc CPU-sec
            % is required for service s, giving this flow proportionality: _sc = d_sc f_c
            % this formula assumes that for every request-response from a client,
            % consumes service rate of the servers proportionally
            % or the service rate at the host level is distributed among the user
            % classes
            gammaV == d .* (ones(o.num_service,1)*X_c')
            X_c >= X_target
            % minimize sum( sum(o.c' * thetaV,2)) +  r1 * pos( X_target - X_c);
            minimize    r(1) * sum(pos( X_target - X_c),1)+...
                r(2) * sum( sum(o.c .* thetaV,2),1);  % r(2) * sum( o.c' * sum(thetaV,2),1);
            cvx_end
            gammaV=gammaV;
            thetaV = round((10^4).*thetaV)/(10^4);
            % recovering theta(o.num_host,o.num_service) from thetaMu
            theta= bsxfun(@rdivide, thetaV, sum(thetaV,1));
            X_c
        end % solve_nfm
        
        % inputs are column wise
        function [f_,rt_, U_ ]=solve_lqm2(o,d_csh,N,Z)
            % here d is  d= rand(num_service,num_classes)
            % the software-hardware demands will be summed up to produce
            % hardware ones
            
            % we add one delay center for each class
            C = size(d_csh,1);
            H=size(d_csh,3);
            K = C+size(d_csh,3);
            S = [diag(Z) permute((sum(d_csh,2)), [2,3,1])];  % think times show up as service demands of delay centers
            V = [eye(C,C) ones(C,H)];  % every class visits his own delay center once
            m = [-1*ones(1,C)  o.cap'];
            [U R Q X] = qncmmva(N, S,V, m );
            % calculating per class measures
            X_= X(1:C, 1:C)*ones(C,1);   % i know V here is 1
            f_= X_;
            R_= sum(R(:,C+1:end) .* V(:,C+1:end), 2);
            rt_ = R_;
            U_= U(:,C+1:end) ;
        end
        
        function [f_,rt_,util_]=solve_lqm(o,theta,N,Z)
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
                    model.nodes=[model.nodes,  OpNode(sprintf('H%d',i),'server', o.cap(i) ,1, 1/o.speed_ratio(i),1)];
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
                                if (d_csh(c,s,h)>0)
                                    scs=model.scenarios;
                                    sc=scs(c);
                                    sc.addCall( OpCall(sprintf('S_c%d',c) , sprintf('S%d_%d',s,h), 1, d_csh(c,s,h), 0)  );
                                end
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
            
            util_=[];
            % The first tool nodes are associated with the classes not the
            % hosts
            for h=1:o.num_host
                util_=[util_ ; model.nodes(h+o.num_classes).cpuUtilization];
            end
        end % solve_lqm
        
        
        % f_sla, T, thetaV_0, cap, r, objective
        function  thetaV= solve_mpc_contention(o,N, Z,X_sla, T, thetaV_0,cap, r, objective)
            X_target = X_sla;
            test=true;
            MAX_ITER= 1;
            thetaV=[];
            X_q=[];
            X_sla 
            LOOKAHEAD_WIND=1;
            for t=1:LOOKAHEAD_WIND
                fprintf('-');
                % disp(sprintf('X_target' ));
                %X_target
                iter=1;
                while true
                    % X_target
                    fprintf('.');
                    [u, thetaV, gammaV, X_nfm,  U_nfm, theta,~] = solve_nfm_over_time...
                        (o, X_target, T, thetaV_0, cap , r, objective);

                    
                    iter=iter+1;
                    if (iter>MAX_ITER) break; end; 
                    X_target_old = X_target;
                    [x_q, r_q] =  solve_lqm(o, theta(:,:,t+1)  , N(:,t+1),  Z);   % this is theta i guess
                    % R_q = [R_q r_q];
                    X_target(:,t+1) = X_target(:,t+1) + (X_sla(:,t+1) - x_q);
                    if ( (X_target(:,t+1)-X_target_old(:,t+1)) < 0.001*X_target(:,t+1)) break; end;
                end % while
                X_q = [X_q x_q];
            end % for
             X_q
            % thetaV
        end % function
        
          function [u_all,thetaV_all, X_nfm_all, U_nfm_all,theta_all,c_infr, c_sla, c_trsh]=...
            iterate_simple_mpc_steps_new(o,X_sla,nsteps,T,r, thetaV_0,objective) 
            u_all=zeros(o.num_host , o.num_service, T);
            thetaV_all=zeros(o.num_host , o.num_service, T+1);
            X_nfm_all=zeros(o.num_classes,T);
            U_nfm_all=zeros(o.num_host,T);
            theta_all=zeros(o.num_host , o.num_service, T+1);
            fprintf('\n running mpc with T=%d',T);
            for t=1:nsteps
                fprintf('.'); 
                X_target_now=X_sla(:,t:t+T-1);      
                [u, thetaV_round, gammaV, X_nfm,  U_nfm, theta,thetaV_orig] = solve_nfm_over_time(...
                    o, X_target_now, T, thetaV_0, o.cap*ones(1,T), r, objective);
               
                
                u_all(:,:,t)=u(:,:,1);
                thetaV_all(: , :, t)=thetaV_round; 
                thetaV_0=thetaV_orig(: , :, 2);
                X_nfm_all(:,t)= x_q; % X_nfm(:,1);
                U_nfm_all(:,t)=U_nfm(:,1);
                theta_all(: , :, t)=  theta(:,:,1);
            end 
            fprintf('\n');
            permute(u_all,[1 3 2]);
            c_infr=sum(o.c' * U_nfm_all);
            c_sla= sum(pos( X_sla(:,2:nsteps)-X_nfm_all(:,2:nsteps)));
            c_trsh= sum(sum(sum(abs(u_all(:,:,:)))));
          end
          
        function [u_all,thetaV_all, X_nfm_all, U_nfm_all,theta_all,c_infr, c_sla, c_trsh]=...
                iterate_simple_mpc_steps(o,X_sla,nsteps,T,r, thetaV_0,objective)  
            
            u_all=zeros(o.num_host , o.num_service, T);
            thetaV_all=zeros(o.num_host , o.num_service, T+1);
            X_nfm_all=zeros(o.num_classes,T);
            U_nfm_all=zeros(o.num_host,T);
            theta_all=zeros(o.num_host , o.num_service, T+1);
            %   thetaV_0 = zeros(o.num_host , o.num_service);
            fprintf('\n running mpc with T=%d',T);
            for t=1:nsteps
                fprintf('.'); % disp(t);
                % this is the perfect information case
                X_target_now=X_sla(:,t:t+T-1); 
                
                %  this is the case where I assume that the workload has
                %  been the same over  the window equal to the current amount
                % X_target_now=repmat(X_sla(:,t), [1 T]);  
                [u, thetaV_round, gammaV, X_nfm,  U_nfm, theta,thetaV_orig] = solve_nfm_over_time(...
                    o, X_target_now, T, thetaV_0, o.cap*ones(1,T), r, objective);
                u_all(:,:,t)=u(:,:,1);
                % this is my state, so i have to follow the system on this
                % can the difference between the the host consumptions and real distribution be state process error
                % and the contention be a observation error?
                thetaV_all(: , :, t)=thetaV_round(:,:,1); % round((10^2).*thetaV_0)/(10^2); % thetaV_0;
                % assume here i solve the real sys to obtain thetaV_0 or x
                 thetaV_0=thetaV_orig(: , :, 2);
                % [i j k]=size(thetaV);   
                % thetaV_0=thetaV(: , :, 2) - rand(i,j).*thetaV(: , :, 2).*0.9; 
                %thetaV_0=thetaV(: , :, 2)-1.1.*thetaV(: , :, 2);
                X_nfm_all(:,t)=X_nfm(:,1);
                U_nfm_all(:,t)=U_nfm(:,1);
                theta_all(: , :, t)=  theta(:,:,1);
                % consumption_cost_tot_all(t)= consumption_cost(1);
            end
            fprintf('\n');
            permute(u_all,[1 3 2]);
            c_infr=sum(o.c' * U_nfm_all);
            % [X_lqm ~]=eveluate_lqm(o, theta_all , N(:,1:nsteps)  ,  Z);
            c_sla= sum(pos( X_sla(:,2:nsteps)-X_nfm_all(:,2:nsteps)));
            c_trsh= sum(sum(sum(abs(u_all(:,:,:)))));
        end
        
        function [X_lqm,r_lqm,util_lqm]=eveluate_lqm(o, theta , N,  Z)
            X_lqm=[]; r_lqm=[]; util_lqm=[];
            for t=1:size(N,2)
                x_q=zeros(size(N,1),1); 
                r_q=zeros(size(N,1),1);
                try
                    [x_q, r_q,util_] =  solve_lqm(o, theta(:,:,t) , N(:,t),  Z);
                catch  err
                    disp(err.message); 
                end
                X_lqm=[X_lqm x_q];
                r_lqm=[r_lqm r_q];
                util_lqm=[util_lqm util_ ];
                %  disp(sprintf('solving lqm for t=%d',t));
                fprintf('.');
            end
        end
        % http://www-stat.wharton.upenn.edu/~stine/stat910/lectures/05_regression.pdf
        %  see
        %  http://www-stat.wharton.upenn.edu/~stine/stat910/lectures/08_intro_arma.pdf
        %  page 15
        function regression(o)
        end
        
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
        
  
    end %methods
end

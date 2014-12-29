
function ex_test()
            num_host=6;
            cap=[8 9 9 8 9 6]; 
            speed_ratio=[1 1.2 0.9 1.1 0.8 1.2];
            num_service= 14;
            num_classes = 2;
          
            Y(1,:)=[1.0000    2.0000    1.0000         0    2.8000    2.5000    1.0000...
                0    5.3000    5.0000    0.5000    1.0000         0    7.9500];
            Y(2,:)=[0         0    2.0000    1.0000         0    7.0000    2.0000    1.0000...
                7.0000   14.0000    1.0000    2.0000    0.5000   10.5000];
            
            D=...
                [0.005, 0.004, 0.002, 0.008, ...
                0.001, 0.002, 0.005, 0.008, ...
                0.006, 0.008, 0.004, 0.005, 0.006,...
                0.005];
            d= Y .* (ones(num_classes,1) * D);
            r_sla=50; 
            r_resource=1;
            r_trsh= 1;
            
            session=ex_create_optimizer(num_host,cap,speed_ratio,num_classes,...
               num_service,r_sla, r_resource, r_trsh,  d); 
           %-------
            
            tmp= workload().get_workload('data/day42_per_min.txt', 1, 24*60);
            tt = ( 0:60:24*60 )';
            int = 1 + floor((1:24*60)./5)';
            mu = accumarray(int,tmp,[],@mean)'; %288
            sd = accumarray(int,tmp,[],@std);
%              N=[0 220 150 200
%                 0 100   170 60];

            nsteps=10; 
            N=[ round(mu(1:nsteps)./5) 
               round(mu(1+nsteps:2*nsteps)./5)];
            Z=[1 1];
            RT_sla=[0.146  0.267];
            T=3;
            N1=[N repmat(N(:,nsteps),[1  T])];
            %X_sla = session.optimizer.calculate_throuput_objective(...
             %   N1, Z, RT_sla, T+nsteps);
            [ X_sla ] = ex_calculate_throuput_objective...
                ( session,N1, Z, RT_sla, T+nsteps )
            t=5; 
           X_target_now=X_sla(:,t:t+T-1); 
            
        thetaV_0 = zeros(num_host , num_service);
        load('theta','thetaV_0')
        for t=1:nsteps
           [theta_,thetaV_0_] = ex_mpc_step(session,X_target_now,T,thetaV_0);  
           thetaV_0=thetaV_0_; 
           %x_q=zeros(size(N,1),1); 
           %r_q=zeros(size(N,1),1);
            try
           [x_q, r_q] =  session.optimizer.solve_lqm( theta_ , N(:,t),  Z)
           catch  err
                    disp(err.message); 
                    disp(t);
           end
           d_csh= permute(  repmat(d, [ 1,1 , num_host] ) , [1,2,3])  .* ...
                  permute( repmat(theta_, [1 , 1, num_classes]), [3,2,1]  );
        end
        
        
%                              % this is going to be more complex for different examples,
%                 % later on, maybe comes from a matrix
%                 % but now its just equal calls from the client to all services
%                 for s=1:o.num_service
%                     for h=1:o.num_host
%                         if (theta(h,s)~=0)
%                             model.services=[model.services,  OpService(sprintf('S%d_%d',s,h), sprintf('T%d_%d',s,h))];
%                             model.containers=[model.containers,   OpContainer(sprintf('T%d_%d',s,h), 1000, sprintf('H%d',h), 'true' )];
%                             for c=1:o.num_classes
%                                 % OpCall(caller, callee, invocations, CPUDemand, DiskDemand)
%                                 if (d_csh(c,s,h)>0)
%                                     scs=model.scenarios;
%                                     sc=scs(c);
%                                     sc.addCall( OpCall(sprintf('S_c%d',c) , sprintf('S%d_%d',s,h), 1, d_csh(c,s,h), 0)  );
%                                 end
%                                 %model.calls=[model.calls...
%                                 %OpCall('ClientS', sprintf('S%d_%d',s,h), 1*theta(h,s), o.d(s,c), 0)];
%                                 % ];
%                             end
%                         end
%                     end
%                 end
                
end

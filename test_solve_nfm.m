 function [thetaV,gammaV]=solve_nfm(d, X_target, num_service, num_classes, num_host, speed_ratio, cap)  
 
           
             r1 = 1000*ones(1,num_classes); 
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
            sum(thetaV,2)<=cap .* speed_ratio             

            sum(thetaV,1)' == sum(gammaV ,2) 
            % For each single user request by class c, a demand of d_sc CPU-sec
            % is required for service s, giving this flow proportionality: _sc = d_sc f_c
            % this formula assumes that for every request-response from a client,
            % consumes service rate of the servers proportionally
            % or the service rate at the host level is distributed among the user
            % classes
             gammaV == d .* (ones(num_service,1)*X_c')  
             %   X_c >= X_target  
            minimize sum( sum(o.c .* thetaV,2)) +  r1 * pos( X_target - X_c); 
            cvx_end 
        end % solve_nfm
        
        function test_solve_nfm1() 
         d =  [0.04 ];  
        cap=[8 8 8 8]';       
         speed_ratio=[1 1 1 1]';             
         num_host=4;  
         num_classes = 1;  
         num_service=1; 
       X_target = 227.2727; 
       solve_nfm(d, X_target, num_service, num_classes, num_host, speed_ratio, cap)
        end 
        
        
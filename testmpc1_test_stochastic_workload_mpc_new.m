classdef testmpc1_test_stochastic_workload_mpc_new
    %UNTITLED8 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
            
        % testmpc1().test_stochastic_workload_mpc_new( 4, 4,50,1, 1, 'abs','opt')
        function [U_nfm,X_nfm,X_sla,relocations,X_lqm,r_lqm,RT_sla,util_lqm]=test_stochastic_workload_mpc_new(...
            obj, nsteps, T,r_sla, r_resource, r_trsh, objective,experiment_kind)
            o=testmpc1();
            testmpc1_service_center_definition.setup_hosts(o);
            testmpc1_service_center_definition.setup_services(o);

            w=workload();
            tmp= w.get_workload('data/day42_per_min.txt', 1, 24*60);
            tt = ( 0:60:24*60 )';
            int = 1 + floor((1:24*60)./5)';
            mu = accumarray(int,tmp,[],@mean)'; %288
            sd = accumarray(int,tmp,[],@std);

%              N=[0 220 150 200
%                 0 100   170 60];
           % N=[ round(mu(1:nsteps)./1)-200 
            %   round(mu(1+nsteps:2*nsteps)./2)-110]; 
            N=[ round(mu(1:nsteps)./5) 
               round(mu(1+nsteps:2*nsteps)./5)];
            Z=[1 1]';
           % RT_sla=[0.146  0.267]';
            RT_sla=[0.17  0.31]';
            RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',10));
             c_orig=rand(o.num_host,1)+10; 
             o.c=c_orig*ones(1, o.num_service);
            % o.c=rand(o.num_host,o.num_service);

            thetaV_0 = zeros(o.num_host , o.num_service);
            r=[r_sla r_resource r_trsh];           
            mpccost = 0; 
            
            if (strcmpi(experiment_kind,'mpc'))
                N1=[N repmat(N(:,nsteps),[1  T])];
                X_sla = calculate_throuput_objective(o, N1, Z, RT_sla, T+nsteps);
                [u,thetaV_round, X_nfm, U_nfm,theta,c_infr, c_sla, c_trsh]=...
                    iterate_simple_mpc_steps(o,X_sla,nsteps,T,r, thetaV_0,objective);
                
            elseif  (strcmpi(experiment_kind,'opt'))
                T=nsteps;
                X_sla = calculate_throuput_objective(o, N, Z, RT_sla, T);
                X_target = X_sla;
                [u, thetaV_round, gammaV, X_nfm,  U_nfm, theta,thetaV_orig] =...
                    solve_nfm_over_time(o, X_target, T, thetaV_0, o.cap*ones(1,T), r, objective);
                %thetaV=round((10^2).*thetaV)/(10^2); 
                %U_host_=permute(sum(thetaV,2),[1 3 2]);
            end
            %-------- this presents the cost of hosts by showing their total number
            % over time
            c_infr=sum(c_orig' *... % sum over time 
            ((0.4*(U_nfm(:,2:nsteps)>0) +...
             0.6*  (U_nfm(:,2:nsteps) ./(o.cap.* o.speed_ratio*ones(1,nsteps-1)))))); 


            %-------- this represents the cost of SLA
            % actually evaluate using lqm 
            [X_lqm, r_lqm, util_lqm]=eveluate_lqm(o, theta , N,  Z);
%            X_nfm=min(X_nfm,X_sla(:,1:nsteps));
            c_sla= sum(sum(pos( X_sla(:,2:nsteps)-X_nfm(:,2:nsteps)),1),2); % X_lqm
%            rt_sla=RT_sla*ones(1,nsteps);
%            c_sla= sum(sum( pos(r_lqm(:,2:nsteps)-rt_sla(:,2:nsteps))  ,1),2);
            
            %-------- this represent the trashing, let's see how we can show it
            % better
            placement = (thetaV_round(:,:,1:nsteps)~=0);    
            pl_delta = placement(:,:,2:nsteps) -  placement(:,:,1:nsteps-1); 
            relocations=permute(sum(sum(pl_delta~=0, 1),2),[1,3,2]);
            %c_trsh= sum(sum(sum(abs(u(:,:,:)))));
            c_trsh= sum(relocations); 

             disp(sprintf('%s_%dsteps_%dT_%d_%d_%d_%s %5.2f, %5.2f, %5.2f',...
                 experiment_kind, nsteps, T,r_sla, r_resource, r_trsh,objective,...
                 c_sla,  c_infr,  c_trsh));
        end
 
              
        %testmpc1().test_test_stochastic_workload_mpc_new()
        function test_test_stochastic_workload_mpc_new(o)
            nsteps=144; 
            T=7;
            r_sla=50; r_resource=10; r_trsh=4; objective='abs'
            [U_nfm,X_nfm,X_sla,relocations,X_lqm,r_lqm,RT_sla,util_lqm]=o.test_stochastic_workload_mpc_new...
                ( nsteps, T,r_sla,r_resource, r_trsh, 'abs','opt')
       
            % subplot(3,1,1); 
            h=figure;             
            plot(sum(U_nfm>0,1));  
            hXLabel=xlabel('Time step'); 
            hYLabel=ylabel('Active Servers');
                      axis([0 nsteps 3 6])
            o.fix_font(hXLabel, hYLabel);
            set(gca, ...
              'Box'         , 'off'     , ...
              'TickDir'     , 'out'     , ...
              'TickLength'  , [.02 .02] , ...
              'XMinorTick'  , 'on'      , ...
              'YMinorTick'  , 'on'      , ...
              'XGrid'       , 'on'      , ...              
              'YGrid'       , 'on'      , ...
              'XColor'      , [.3 .3 .3], ...
              'YColor'      , [.3 .3 .3], ...
              'XTick'       , 0:20:144, ...
              'LineWidth'   , 2         ,...
              'YTick'       , 0:1:7, ...
              'LineWidth'   , 2         );
            try
            UtilityLib.print_figure(h,9,7,sprintf('figure/servers_%dsteps_%dT_%d_%d_%d_%s',...
                            nsteps, T,r_sla, r_resource, r_trsh,objective)); 
            catch me;  me;  end
               
             %subplot(3,1,2); 
            h=figure;
            %plot(pos( X_sla(:,2:T)-X_nfm_opt(:,2:T)));
            %rt_sla=(N./X_sla)-Z; rt=(N./X_nfm_opt)-Z;
            %plot(2:T, rt_sla(:,2:T),'-.b',2:T,rt(:,2:T),'-r');
            %[f,x] = ecdf(X_nfm(1,2:T)-X_sla(1,2:T));  plot(x,f,'r'); hold on;
            %[f,x] = ecdf(X_nfm(2,2:T)-X_sla(2,2:T));  plot(x,f,'b');
            %p=cdfplot(X_nfm(1,2:T)-X_sla(1,2:T)); set(h,'Color','r');hold on;
            %p=cdfplot(X_nfm(2,2:T)-X_sla(2,2:T)); set(h,'Color','b');
            [f,x] = ecdf( r_lqm(1,2:nsteps) - repmat(RT_sla(1),[1 nsteps-1])); 
            plot(x,f,'r'); hold on;
            [f,x] = ecdf( r_lqm(2,2:nsteps) - repmat(RT_sla(2),[1 nsteps-1])); 
            plot(x,f,'b');
            hLegend = legend('$RT_{c1}-RT^{SLA}_{c1}$','$RT_{c2}-RT^{SLA}_{c2}$');
            set([hLegend, gca],'FontSize', 18 );          
            set(hLegend,'Interpreter','Latex');
            hXLabel=xlabel('Response Time Mismatch'); 
            hYLabel=ylabel('Cumulative Probablity');  
            try
            UtilityLib.print_figure(h,9,7,sprintf('figure/rtcdf_%dsteps_%dT_%d_%d_%d_%s',...
                nsteps, T,r_sla, r_resource, r_trsh,objective)); 
            catch me;  me;  end
%             plot(2:nsteps, X_sla(1,2:nsteps),'-r',...
%                 2:nsteps, X_nfm(1,2:nsteps),'--r',...
%                 2:nsteps, X_sla(2,2:nsteps),'-b',...
%                 2:nsteps, X_nfm(2,2:nsteps),'-.b');  
%             hLegend = legend('$X^{SLA}_{c1}$',...
%                            '$X_{c1}$',...
%                            '$X^{SLA}_{c2}$',...
%                            '$X_{c2}$');
             h=figure; 
             plot(2:nsteps, repmat(RT_sla(1),[1 nsteps-1]),'-r',...
              2:nsteps, r_lqm(1,2:nsteps),'--r',...
              2:nsteps, repmat(RT_sla(2),[1 nsteps-1]),'-b',...
              2:nsteps, r_lqm(2,2:nsteps),'-.b');      
            hLegend = legend('$RT^{SLA}_{c1}$',...
                           '$RT_{c1}$',...
                           '$RT^{SLA}_{c2}$',...
                           '$RT_{c2}$');
            set([hLegend, gca],'FontSize', 20 );          
            set(hLegend,'Interpreter','Latex');
            hXLabel=xlabel('Time step'); 
            hYLabel=ylabel('Response Time');  
           try
            UtilityLib.print_figure(h,9,7,sprintf('figure/rt1_%dsteps_%dT_%d_%d_%d_%s',...
                nsteps, T,r_sla, r_resource, r_trsh,objective)); 
           catch me;  me;  end
           
            % subplot(3,1,3); 
            h=figure; 
            %plot(relocations);  % over hosts,then services
            [f,x] = ecdf(relocations);  plot(x,f); 
            axis([0 20 0 1]);
            hXLabel=xlabel('Relocation per Step'); 
            hYLabel=ylabel('Cumulative Probablity'); 
            o.fix_font(hXLabel, hYLabel);
              set(gca, ...
              'Box'         , 'off'     , ...
              'TickDir'     , 'out'     , ...
              'TickLength'  , [.02 .02] , ...
              'XMinorTick'  , 'on'      , ...
              'YMinorTick'  , 'on'      , ...
              'XGrid'       , 'on'      , ...              
              'YGrid'       , 'on'      , ...
              'XColor'      , [.3 .3 .3], ...
              'YColor'      , [.3 .3 .3], ...
              'XTick'       , 0:10:100, ...
              'LineWidth'   , 2         ,...
              'YTick'       , 0:0.1:1, ...
              'LineWidth'   , 2         );
          try
           UtilityLib.print_figure(h,9,7,sprintf('figure/relocations_%dsteps_%dT_%d_%d_%d_%s',...
                nsteps, T,r_sla, r_resource, r_trsh,objective));   
          catch me;  me;  end
          
           h=figure; 
           plot(util_lqm');
            hXLabel=xlabel('Time Step'); 
            hYLabel=ylabel('Utilization'); 
            o.fix_font(hXLabel, hYLabel);
           set([hLegend, gca],'FontSize', 20 );          
           set(hLegend,'Interpreter','Latex');
           legend('H1','H2','H3','H4','H5','H6'); 
           try
           UtilityLib.print_figure(h,9,7,sprintf('figure/util_%dsteps_%dT_%d_%d_%d_%s',...
                nsteps, T,r_sla, r_resource, r_trsh,objective));
           catch me;  me;  end
        end
  
             %  testmpc1().test_stochastic_workload_new( 100, 5, 4000, 'square')
        %  testmpc1().test_stochastic_workload_new( 144, 1000,50,1, 1, 'abs')
        function test_stochastic_workload_new(...
                o, nsteps, T,r_sla, r_resource, r_trsh, objective)
             o.setup_hosts();
            o.setup_services();
%             N=[];
%             N=[0 220 150 200
%                  0 100   170 60];
            w=workload();
            tmp= w.get_workload('data/day42_per_min.txt', 1, 24*60);
            tt = ( 0:60:24*60 )';
            int = 1 + floor((1:24*60)./5)';
            mu = accumarray(int,tmp,[],@mean)'; %288
            sd = accumarray(int,tmp,[],@std);
            
            N=[ round(mu(1:nsteps)./1)-200 
                round(mu(1+nsteps:2*nsteps)./2)-110];           
            Z=[1 1]';
            RT_sla=[0.146  0.267]';
            %h=figure;   
            plot(1:nsteps, N(1,:) ,'-r'  , 'LineWidth',2); hold on;
            plot(1:nsteps, N(2,:) ,'-b'  , 'LineWidth',2);  
            hLegend = legend('$N_{c1}$','$N_{c2}$');
            set(hLegend,'Interpreter','Latex');
            set([hLegend, gca],'FontSize', 24 );
            % title('Workload');
            hXLabel = xlabel('Time step'); 
            hYLabel = ylabel('Number of Users (N)');  
            o.fix_font(hXLabel, hYLabel);
            
            set(gca, ...
              'Box'         , 'off'     , ...
              'TickDir'     , 'out'     , ...
              'TickLength'  , [.02 .02] , ...
              'XMinorTick'  , 'on'      , ...
              'YMinorTick'  , 'on'      , ...
              'XGrid'       , 'on'      , ...              
              'YGrid'       , 'on'      , ...
              'XColor'      , [.3 .3 .3], ...
              'YColor'      , [.3 .3 .3], ...
              'XTick'       , 0:20:144, ...
              'LineWidth'   , 2         ,...
              'YTick'       , 0:100:400, ...
              'LineWidth'   , 2         );
            %UtilityLib.print_figure(h,9,7,sprintf('figure/workload_%dsteps_%dT_%d_%d_%d_%s',...
            %    nsteps, T,r_sla, r_resource, r_trsh,objective));

            RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',10));
             c_orig=rand(o.num_host,1)+10; 
             o.c=c_orig*ones(1, o.num_service);
            % o.c=rand(o.num_host,o.num_service);

            thetaV_0 = zeros(o.num_host , o.num_service);
            r=[r_sla r_resource r_trsh];           
            mpccost = 0; 


%             N1=[N repmat(N(:,nsteps),[1  T])];
%             X_sla = calculate_throuput_objective(o, N1, Z, RT_sla, T+nsteps);
%             [u,thetaV, X_nfm, U_nfm,theta,c_infr, c_sla, c_trsh]=...
%                 iterate_simple_mpc_steps(o,X_sla,nsteps,T,r, thetaV_0,objective);

            T=nsteps;
            X_sla = calculate_throuput_objective(o, N, Z, RT_sla, T);
            X_target = X_sla;
            [u, thetaV, gammaV, X_nfm,  U_nfm, theta] =...
                solve_nfm_over_time(o, X_target, T, thetaV_0, o.cap*ones(1,T), r, objective);

            % this presents the cost of hosts by showing their total number
            % over time
            
            c_infr=sum(c_orig' *... % sum over time 
            ((0.4*(U_nfm(:,2:T+1)>0) +...
             0.6*  (U_nfm(:,2:T+1) ./(o.cap.* o.speed_ratio*ones(1,T)))))); 
            % o.c' *(U_nfm>0);
            % subplot(3,1,1); 
            %h=figure
            plot(sum(U_nfm>0,1));  
            hXLabel=xlabel('Time step'); 
            hYLabel=ylabel('Active Servers');
            axis([0 nsteps 0 6])
            o.fix_font(hXLabel, hYLabel);
            set(gca, ...
              'Box'         , 'off'     , ...
              'TickDir'     , 'out'     , ...
              'TickLength'  , [.02 .02] , ...
              'XMinorTick'  , 'on'      , ...
              'YMinorTick'  , 'on'      , ...
              'XGrid'       , 'on'      , ...              
              'YGrid'       , 'on'      , ...
              'XColor'      , [.3 .3 .3], ...
              'YColor'      , [.3 .3 .3], ...
              'XTick'       , 0:20:144, ...
              'LineWidth'   , 2         ,...
              'YTick'       , 0:1:7, ...
              'LineWidth'   , 2         );
            %UtilityLib.print_figure(h,9,7,sprintf('figure/active_servers_%dsteps_%dT_%d_%d_%d_%s',...
            %    nsteps, T,r_sla, r_resource, r_trsh,objective));

            % here am planning to have one response time
            % cumulative distribution
            % function(CDF) for each class on the same plot
            % representing with the SLA as  horizontal line
            % the CDF after the horizontal line should not grow much
            % 
            % we do not evaluate for the steps we do not control
            % c_sla_opt= sum(pos( X_sla(:,2:T)-X_nfm(:,2:T)));
            % [X_lqm ~]=eveluate_lqm(o, theta , N,  Z);
            %c_sla_opt= sum(pos( X_sla(:,2:T)-X_nfm_opt(:,2:T)));
            c_sla= sum(sum(pos( X_sla(:,2:T)-X_nfm(:,2:T)),1),2);
            %subplot(3,1,2); 
            %h=figure
            %plot(pos( X_sla(:,2:T)-X_nfm_opt(:,2:T)));
            %rt_sla=(N./X_sla)-Z; rt=(N./X_nfm_opt)-Z;
            %plot(2:T, rt_sla(:,2:T),'-.b',2:T,rt(:,2:T),'-r');
            [f,x] = ecdf(X_nfm(1,2:T)-X_sla(1,2:T));  plot(x,f,'r'); hold on;
            [f,x] = ecdf(X_nfm(2,2:T)-X_sla(2,2:T));  plot(x,f,'b');
            %p=cdfplot(X_nfm(1,2:T)-X_sla(1,2:T)); set(h,'Color','r');hold on;
            %p=cdfplot(X_nfm(2,2:T)-X_sla(2,2:T)); set(h,'Color','b');
            hLegend = legend('$X_{c1}-X^{SLA}_{c1}$','$X_{c2}-X^{SLA}_{c2}$');
            set([hLegend, gca],'FontSize', 18 );          
            set(hLegend,'Interpreter','Latex');
            hXLabel=xlabel('Throughput Mismatch'); 
            hYLabel=ylabel('Cumulative Probablity');  
%             plot(2:T, X_sla(1,2:T),'-r',...
%                 2:T, X_nfm(1,2:T),'--r',...
%                 2:T, X_sla(2,2:T),'-b',...
%                 2:T, X_nfm(2,2:T),'-.b');     
%             hLegend = legend('$X^{SLA}_{c1}$',...
%                            '$X_{c1}$',...
%                            '$X^{SLA}_{c2}$',...
%                            '$X_{c2}$');
%             set([hLegend, gca],'FontSize', 20 );          
%             set(hLegend,'Interpreter','Latex');
%             hXLabel=xlabel('Time step'); 
%             hYLabel=ylabel('Throughput');  
            o.fix_font(hXLabel, hYLabel);
              set(gca, ...
              'Box'         , 'off'     , ...
              'TickDir'     , 'out'     , ...
              'TickLength'  , [.02 .02] , ...
              'XMinorTick'  , 'on'      , ...
              'YMinorTick'  , 'on'      , ...
              'XGrid'       , 'on'      , ...              
              'YGrid'       , 'on'      , ...
              'XColor'      , [.3 .3 .3], ...
              'YColor'      , [.3 .3 .3], ...
              'XTick'       , 0:20:144, ...
              'LineWidth'   , 2         ,...
              'YTick'       , 0:0.1:1, ...
              'LineWidth'   , 2         );
            %UtilityLib.print_figure(h,9,7,sprintf('figure/throughcdf_%dsteps_%dT_%d_%d_%d_%s',...
            %    nsteps, T,r_sla, r_resource, r_trsh,objective));

           
            % this represent the trashing, let's see how we can show it
            % better
            placement = (thetaV(:,:,1:T)~=0);    
            pl_delta = placement(:,:,2:T) -  placement(:,:,1:T-1); 
            relocations=permute(sum(sum(pl_delta~=0, 1),2),[1,3,2]);
            %c_trsh= sum(sum(sum(abs(u(:,:,:)))));
            c_trsh= sum(relocations); 
            % subplot(3,1,3); 
            %h=figure 
            %plot(permute(sum(sum(pl_delta~=0, 1),2),[1,3,2]));  % over hosts,then services
            [f,x] = ecdf(relocations);  plot(x,f); 
%             edges = linspace(0,30,100);
%             [count] = histc(relocations , edges )';
%             bar( edges, count, 'histc' );
%             hline = findobj(gca,'Type','line'); delete(hline)
%             hpatch = findobj(gca,'Type','patch');
%             set(hpatch,'FaceColor',[0.91,0.91,0.91])
            hXLabel=xlabel('Relocation per Step'); 
            hYLabel=ylabel('Cumulative Probablity'); 
            o.fix_font(hXLabel, hYLabel);
              set(gca, ...
              'Box'         , 'off'     , ...
              'TickDir'     , 'out'     , ...
              'TickLength'  , [.02 .02] , ...
              'XMinorTick'  , 'on'      , ...
              'YMinorTick'  , 'on'      , ...
              'XGrid'       , 'on'      , ...              
              'YGrid'       , 'on'      , ...
              'XColor'      , [.3 .3 .3], ...
              'YColor'      , [.3 .3 .3], ...
              'XTick'       , 0:10:100, ...
              'LineWidth'   , 2         ,...
              'YTick'       , 0:0.1:1, ...
              'LineWidth'   , 2         );
            %UtilityLib.print_figure(h,9,7,sprintf('figure/relocations_%dsteps_%dT_%d_%d_%d_%s',...
            %    nsteps, T,r_sla, r_resource, r_trsh,objective));

             disp(sprintf('%dsteps_%dT_%d_%d_%d_%s %5.2f, %5.2f, %5.2f',...
                 nsteps, T,r_sla, r_resource, r_trsh,objective,...
                 c_sla,  c_infr,  c_trsh));
%             optcost =  r(1)*c_sla_opt + r(2) * c_infr_opt + r(3)*c_trsh_opt
% 
%             figure; 
%             edges = linspace(0,0.6,100);
%             [count] = histc( permute(sum(abs(u_opt(:,:,:))), [1 3 2]), edges )';
%             bar( edges, count, 'histc' );
%             hline = findobj(gca,'Type','line'); delete(hline)
%             hpatch = findobj(gca,'Type','patch');
%             set(hpatch,'FaceColor',[0.91,0.91,0.91])
        end
    
          
      
% this is the part that I actually  do not use the li's example

        % [optcost_tr, costmpc, optcost] = testmpc1().test_stochastic_workload( 100, 5, 4000, 'square')
        function [optcost_trIg, mpccost, optcost] = test_stochastic_workload(...
                o, nsteps, T,r_trsh_cost, objective)
            o.num_host=8;
            o.num_classes = 1;
            o.num_service=1;
            o.cap=[8 8 8 8 8 8 8 8]';
            o.speed_ratio=[1 1 1 1 1 1 1 1]';
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
            D =  [0.14 ];
            o.d= Y .* (ones(o.num_classes,1) * D);
            Z=[1]';
            RT_sla=[0.146]';
             r=[50 1 r_trsh_cost];
             
            %   % mpc d
            % nsteps = 10;
            %  T=5;
            thetaV_0 =   [0   8.0000   0.0360  0 0 0 0 0]';
           
            w=workload();
            tmp= w.get_workload('data/day42_per_min.txt', 1, 24*60)';
            % N_orig=tmp(1040:1040+nsteps-1+10);  
            N_orig=tmp(1:1+nsteps-1+10); 
            a = 1;
            % b = [1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8];
            b = ones(1,8)/8; 
            N = filter(b,a,N_orig);  
            N=N(11:end);
            N_orig=N_orig(11:end); 
%              h=figure;
%             plot(1:nsteps, N_orig,'--b'); hold on; 
%             plot(1:nsteps, N ,'-r'  , 'LineWidth',2);  
%             hleg1 = legend('N','N_{smooth}' );
%             title('Workload');
%             xlabel('Time step'); ylabel('Number of Users (N)');  
%             axis([0 nsteps min(N_orig)-20 max(N_orig)+20])
%             UtilityLib.print_figure(h,9,7,sprintf('figure/workload'));
            
           % N=200+[1:nsteps]+50*sin([1:5:5*(nsteps)]/20);
           
            N1=[N ones(1,T)*N(nsteps)];
            % plot(N); hold on;
            X_sla = calculate_throuput_objective(o, N1, Z, RT_sla, T+nsteps);
           
%             [u_mpc,thetaV_mpc, X_nfm_mpc, U_nfm_mpc,theta_mpc,c_infr, c_sla, c_trsh]=...
%                 iterate_simple_mpc_steps(o,X_sla,nsteps,T,r, thetaV_0,objective);
%             disp(sprintf('\n c_mpc_sla=%5.4f,  c_mpc_infr=%5.4f, c_mpc_trsh=%5.4f', c_sla,  c_infr, c_trsh));
%             mpccost =  r(1)*c_sla + r(2) * c_infr+ r(3)*c_trsh
%             permute(sum(abs(u_mpc(:,:,:))), [1 3 2]) ;
            mpccost = 0; 
             
            % this one igores the cost of trashing totally,
            % so it is off
            T_orig=T;
            T=nsteps;
            % N=200+[1:T]+50*sin([1:5:5*T]/10);
            X_sla = calculate_throuput_objective(o, N, Z, RT_sla, T);
            % thetaV_0 = zeros(o.num_host , o.num_service);
            % thetaV_0 =   [0   8.0000   0.0360  0]';
            X_target = X_sla;
            r_trsh=[50 1 0.01];
%             [u_opt_trIg, thetaV_trIg, gammaV_trIg, X_nfm_trIg,  U_nfm_trIg, theta_trIg] =...
%                 solve_nfm_over_time(o, X_target, T, thetaV_0, o.cap*ones(1,T), r_trsh, objective);
%             c_infr_opt_trIg=sum(o.c' * U_nfm_trIg);
%             % we do not evaluate for the steps we do not control
%             % [X_lqm ~]=eveluate_lqm(o, theta , N,  Z);
%             c_sla_opt_trIg= sum(pos( X_sla(:,2:T)-X_nfm_trIg(:,2:T)));
%             c_trsh_opt_trIg= sum(sum(sum(abs(u_opt_trIg(:,:,:)))));
%             disp(sprintf('\n c_sla_opt_trIg=%5.4f,  c_infr_opt_trIg=%5.4f,  c_trsh_opt_trIg=%5.4f,', ...
%                 c_sla_opt_trIg,  c_infr_opt_trIg,  c_trsh_opt_trIg));
%             optcost_trIg =  r(1)*c_sla_opt_trIg + r(2) * c_infr_opt_trIg + r(3)*c_trsh_opt_trIg
%             permute(sum(abs(u_opt_trIg(:,:,:))), [1 3 2]);
            optcost_trIg=0; 
            
            % this is the optimal for the simplistic approach evaluated
            % with lqm
            [u_opt, thetaV, gammaV, X_nfm_opt,  U_nfm_opt, theta_opt] =...
                solve_nfm_over_time(o, X_target, T, thetaV_0, o.cap*ones(1,T), r, objective);
            %permute(u_opt,[1 3 2]);
            figure;
            % this presents the cost of hosts by showing their total number
            % over time
            %c_infr_opt=sum(o.c' * U_nfm_opt);
            c_infr_opt=o.c' *(U_nfm_opt>0);
            subplot(3,1,1); plot(sum(U_nfm_opt>0,1));
            
            % here am planning to have one response time
            % cumulative distribution
            % function(CDF) for each class on the same plot
            % representing with the SLA as  horizontal line
            % the CDF after the horizontal line should not grow much
            % 
            % we do not evaluate for the steps we do not control
            % c_sla_opt= sum(pos( X_sla(:,2:T)-X_nfm(:,2:T)));
            % [X_lqm ~]=eveluate_lqm(o, theta , N,  Z);
            %c_sla_opt= sum(pos( X_sla(:,2:T)-X_nfm_opt(:,2:T)));
            c_sla_opt= sum(pos( X_sla(:,2:T)-X_nfm_opt(:,2:T)));
            subplot(3,1,2); 
            %plot(pos( X_sla(:,2:T)-X_nfm_opt(:,2:T)));
            %rt_sla=(N./X_sla)-Z; rt=(N./X_nfm_opt)-Z;
            %plot(2:T, rt_sla(:,2:T),'-.b',2:T,rt(:,2:T),'-r');

%             edges = linspace(-10,10,100);
%             [count] = histc( X_nfm_opt(:,2:T)-X_sla(:,2:T), edges )';
%             bar( edges, count, 'histc' );
%             hline = findobj(gca,'Type','line'); delete(hline)
%             hpatch = findobj(gca,'Type','patch');
%             set(hpatch,'FaceColor',[0.91,0.91,0.91])
             cdfplot(X_nfm_opt(:,2:T)-X_sla(:,2:T))
            
            % this represent the trashing, let's see how we can show it
            % better
            %c_trsh_opt= sum(sum(sum(abs(u_opt(:,:,:)))));
            c_trsh_opt= sum(sum(sum(abs(u_opt(:,:,:)))));
            placement = (thetaV(:,:,1:T)~=0);    
            pl_delta = placement(:,:,2:T) -  placement(:,:,1:T-1);  
            subplot(3,1,3);
             plot(permute(sum(sum(pl_delta~=0, 1),2),[1,3,2]));  % over hosts,then services
%             edges = linspace(0,10,100);
%             [count] = histc( permute(sum(sum(pl_delta~=0, 1),2),[1,3,2]), edges )';
%             bar( edges, count, 'histc' );
%             hline = findobj(gca,'Type','line'); delete(hline)
%             hpatch = findobj(gca,'Type','patch');
%             set(hpatch,'FaceColor',[0.91,0.91,0.91])

            
            disp(sprintf('\n c_sla_opt=%5.4f, c_infr_opt=%5.4f, c_trsh_opt=%5.4f',...
                c_sla_opt,  c_infr_opt,  c_trsh_opt));
            optcost =  r(1)*c_sla_opt + r(2) * c_infr_opt + r(3)*c_trsh_opt
            
%             % compare cost of infrastructure
%             %  plot(1:T, o.c' *U_nfm_tr,'b-', 1:T, o.c' *U_nfm_notr ,'r--' );
%             h=figure;
%             plot(1:T, o.c' *U_nfm_trIg,'-.b',...
%                     1:T, o.c' *U_nfm_opt ,'-r',... 
%                     1:T, o.c' *U_nfm_mpc, '--g'); 
%             hleg1 = legend('one step','optimal' ,'mpc');
%             title('Cost of infrastructure');
%             xlabel('Time step'); ylabel('Infrastructure Cost');
%             UtilityLib.print_figure(h,9,7,sprintf('figure/infrastructure_cost_nsteps%d_T%d_rTRSH%d%s', nsteps, T_orig,r_trsh_cost,objective)); 

%             % compare cost of trashing
%             h=figure;
%             plot(1:T, permute(sum(abs(u_opt_trIg(:,:,:))), [1 3 2]),  '-.b',...
%                 1:T,  permute(sum(abs(u_opt(:,:,:))), [1 3 2]),'-r',...
%                 1:T, permute(sum(abs(u_mpc(:,:,:))), [1 3 2])  , '--g');
%             hleg1 = legend('one step','optimal' ,'mpc');
%             title('Cost of trashing');
%             xlabel('Time step'); ylabel('Trashing Cost');
%             UtilityLib.print_figure(h,9,7,sprintf('figure/trashing_cost_nsteps%d_T%d_rTRSH%d%s', nsteps, T_orig,r_trsh_cost,objective));
        
%            % histogram plots
%             figure(2), clf
%             subplot(3,1,1); 
%             edges = linspace(0,0.6,100);
%             [count] = histc( permute(sum(abs(u_opt_trIg(:,:,:))), [1 3 2]), edges )';
%             bar( edges, count, 'histc' );
%             hline = findobj(gca,'Type','line'); delete(hline)
%             hpatch = findobj(gca,'Type','patch');
%             set(hpatch,'FaceColor',[0.91,0.91,0.91])

%            subplot(3,1,2); 
figure; 
            edges = linspace(0,0.6,100);
            [count] = histc( permute(sum(abs(u_opt(:,:,:))), [1 3 2]), edges )';
            bar( edges, count, 'histc' );
            hline = findobj(gca,'Type','line'); delete(hline)
            hpatch = findobj(gca,'Type','patch');
            set(hpatch,'FaceColor',[0.91,0.91,0.91])

            
%             subplot(3,1,3); 
%             edges = linspace(0,0.05,100);
%             [count] = histc( permute(sum(abs(u_mpc(:,:,:))), [1 3 2]), edges )';
%             bar( edges, count, 'histc' );
%             hline = findobj(gca,'Type','line'); delete(hline)
%             hpatch = findobj(gca,'Type','patch');
%             set(hpatch,'FaceColor',[0.91,0.91,0.91])
%             % axis([-2 2 0 150]);
%             %print -depsc online_reg_dist
        end
   
            % testmpc1().generate_two_CDF()
        function generate_two_CDF(o)
            nsteps=144;
            T=4;
            r_sla=50; 
            objective='abs';
            
%             r_resource=10; r_trsh=4; 
%             [U_nfm1,X_nfm1,X_sla1,relocations1]=testmpc1().test_stochastic_workload_mpc_new...
%                 ( nsteps, T,r_sla,r_resource, r_trsh, 'abs','opt')
% 
%             r_resource=1; r_trsh=40; 
%             [U_nfm2,X_nfm2,X_sla2,relocations2]=testmpc1().test_stochastic_workload_mpc_new...
%                 ( nsteps, T,r_sla,r_resource, r_trsh, 'abs','opt')
            
            
            h=figure; 
            %plot(relocations);  % over hosts,then services
            [f,x] = ecdf(relocations1);  plot(x,f,'r'); hold on;
            [f,x] = ecdf(relocations2);  plot(x,f,'--b'); 
            axis([0 20 0 1]);
            hXLabel=xlabel('Relocation per Step'); 
            hYLabel=ylabel('Cumulative Probablity'); 
            o.fix_font(hXLabel, hYLabel);
              set(gca, ...
              'Box'         , 'off'     , ...
              'TickDir'     , 'out'     , ...
              'TickLength'  , [.02 .02] , ...
              'XMinorTick'  , 'on'      , ...
              'YMinorTick'  , 'on'      , ...
              'XGrid'       , 'on'      , ...              
              'YGrid'       , 'on'      , ...
              'XColor'      , [.3 .3 .3], ...
              'YColor'      , [.3 .3 .3], ...
              'XTick'       , 0:10:100, ...
              'LineWidth'   , 2         ,...
              'YTick'       , 0:0.1:1, ...
              'LineWidth'   , 2         );
           UtilityLib.print_figure(h,9,7,sprintf('figure/relocations_CDF_compare',...
                nsteps, T,r_sla, r_resource, r_trsh,objective));

            h=figure
            plot(sum(U_nfm1>0,1),'r');  hold on;
            plot(sum(U_nfm2>0,1),'--b'); 
            hXLabel=xlabel('Time step'); 
            hYLabel=ylabel('Active Servers');
            axis([0 nsteps 0 6])
            o.fix_font(hXLabel, hYLabel);
            set(gca, ...
              'Box'         , 'off'     , ...
              'TickDir'     , 'out'     , ...
              'TickLength'  , [.02 .02] , ...
              'XMinorTick'  , 'on'      , ...
              'YMinorTick'  , 'on'      , ...
              'XGrid'       , 'on'      , ...              
              'YGrid'       , 'on'      , ...
              'XColor'      , [.3 .3 .3], ...
              'YColor'      , [.3 .3 .3], ...
              'XTick'       , 0:20:144, ...
              'LineWidth'   , 2         ,...
              'YTick'       , 0:1:7, ...
              'LineWidth'   , 2         );
            UtilityLib.print_figure(h,9,7,sprintf('figure/active_servers_CDF_compare',...
                nsteps, T,r_sla, r_resource, r_trsh,objective));
        end
        
        % testmpc1().test_mpc_horizon()
       function test_mpc_horizon(o)
            for i=1:10
                 testmpc1().test_stochastic_workload_mpc_new( 144, i ,50,1, 5, 'abs','mpc');
            end
       end
 
       function draw_mpc_horizon(o)
            r_sla=50;
            r_resource=1;
            r_trashing=1;
            objective='abs';
            data=[23537.72, 50.32,  0.00
                0.95, 8878.70, 24.00
                0.95, 8340.34, 24.00
                0.95, 7602.69, 30.00
                0.95, 6984.34, 59.00
                0.95, 6938.55, 56.00
                0.95, 6962.87, 52.00
                0.95, 6981.97, 60.00
                0.95, 6997.69, 52.00
                0.95, 6998.63, 50.00]; 
            data(:,2)= data(:,2)./10;
            f=figure;   
            h=bar(data);
            hLegend=legend('c_{SLA}',  'c_{infr}',  'c_{trsh}');
            x=get(h,'Xdata');
            y=get(h,'Ydata');
            for xc=1:size(data,1)
                text(xc,max(data(xc,:))+20,[num2str(xc)]);
            end
            UtilityLib.print_figure(f,9,7,sprintf('figure/horizon_%d_%d_%d_%s',...
                r_sla, r_resource,r_trashing, objective));
       end
       
       function fix_font(o,hXLabel, hYLabel )
            set( gca                       , ...
                'FontName'   , 'Times New Roman' );
            set([hXLabel, hYLabel], ...
                'FontName'   , 'Times New Roman');
            set([hXLabel, hYLabel]  , ...
                'FontSize'   , 24          );
         
        end
        
     
    end
    
end


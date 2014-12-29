classdef testmpc1_test_trade_off
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
          function test_trade_off(o,T)
            % ratios=[1 5  20 50 100 200 400 800 1600 3200];
            ratios=[0]; 
            for i=1:10
                % testmpc1().test_stochastic_workload_new( 144, 5, 1000,1, ratios(i), 'abs');
                testmpc1().test_stochastic_workload_mpc_new( 144, 1 ,50,1, ratios(i), 'abs','mpc');
            end
        end
        
        function test_trade_off_and_horizon(o)
            ratios=[1 5  20 50 100 200 400 800 1600 3200];
            for horizon=3:10
                for cost_index=1:10
                    try
                        % testmpc1().test_stochastic_workload_new( 144, 5, 1000,1, ratios(i), 'abs');
                        testmpc1().test_stochastic_workload_mpc_new( 144, horizon ,50,1, ratios(cost_index), 'abs','mpc');
                    catch ME
                        ME
                    end
                end
            end
        end
        
        function draw_trade_off_curves(o)
             %50,1, 
              r_sla=50;
              r_resource=1;
              objective='abs';
             ratios=[1 5  20 50 100 200 400 800 1600 3200];
             data_rSLA50=[...
             0.95, 5440.00, 756.00
              0.95, 5676.41, 388.00
              0.95, 6037.13, 236.00
              11.70, 6485.29, 130.00
              27.75, 6572.64, 130.00
              118.73, 6910.66, 88.00
              193.45, 6914.61, 70.00
              380.00, 6895.28, 70.00
              585.15, 6849.88, 65.00];
          data_rSLA1000=[...
              0.95, 5439.68, 756.00
              0.95, 5440.09, 756.00
              0.95, 5676.29, 388.00
              0.95, 6036.76, 236.00
              0.95, 6487.03, 134.00
              0.95, 6623.52, 130.00
              0.95, 7102.56, 92.00
              0.95, 7371.85, 92.00
              0.95, 8314.25, 84.00
              0.95, 8323.31, 84.00];
          %50 1 (5:14)
          data_50_1_5to14_MPC2=[...
              0.95, 5685.81, 788.00
              0.95, 5629.86, 792.00
              0.95, 5647.02, 774.00
              0.95, 5625.50, 782.00
              0.95, 5660.01, 808.00
              0.95, 5638.53, 828.00
              0.95, 8301.08, 95.00
              0.95, 8301.83, 95.00
              0.95, 8301.60, 95.00
              0.95, 8301.64, 95.00];  
          data=data_50_1_5to14_MPC2;
          data(:,2)= data(:,2)./10;
            f=figure; 
            h=bar(data);
            % hLegend=legend('c_{SLA}',  'c_{infr}',  'c_{trsh}');
            hXLabel=xlabel('trashing cost coefficient index'); 
            hYLabel=ylabel('horizon'); 
            o.fix_font(hXLabel, hYLabel);
            x=get(h,'Xdata');
            y=get(h,'Ydata');
            for xc=1:size(data,1)
                %text(x{2}(xc),y{2}(xc),['y=' num2str(ratios(xc))]);
                text(xc,max(data(xc,:))+20,[num2str(ratios(xc))]);
            end
            
            set(gca, ...
              'Box'         , 'off'     , ...
              'TickDir'     , 'out'     , ...
              'TickLength'  , [.02 .02] , ...
              'XMinorTick'  , 'on'      , ...
              'YMinorTick'  , 'on'      , ...
              'XGrid'       , 'on'      , ...              
              'YGrid'       , 'on'      , ...
              'XColor'      , [.3 .3 .3], ...
              'YColor'      , [.3 .3 .3], ... %'XTick'       , 0:20:9, ...
              'LineWidth'   , 2         ,...
              'YTick'       , 0:100:1000, ...
              'LineWidth'   , 2         );
          
            %UtilityLib.print_figure(f,9,7,sprintf('figure/trade_off_%d_%d_%s',...
                %r_sla, r_resource,objective));
        end
        
        function draw_trade_off_and_horizon_surface(o)
              cost_elements=testmpc1().data_trade_off_and_horizon();
              % cost_elements=cost_elements(1:8,:,:); 
              cost_elements=permute( cost_elements,[2,3,1]); 
              cost_elements(:,2,:)= cost_elements(:,2,:)./10; 
              cost_elements=cost_elements(2:end,:,:);
               % bar(cost_elements(:,3,1))
              for i= 9:-1:1
                tt=permute( cost_elements(i,:,:),[3 2 1])
                scatter3(tt(2:10,1),tt(2:10,2),...
                    tt(2:10,3),40*ones(9,1),repmat([i/10 0 0],[9 1])); 
                axis([0 900 0 900 0 900]); 
                 hXLabel=xlabel('SLA cost'); 
                 hYLabel=ylabel('resource cost'); 
                 hZLabel=zlabel('trashing cost'); 
               hold on;
              end
        end
        
        function draw_trade_off_and_horizon(o)
            cost_elements=testmpc1().data_trade_off_and_horizon();
%             cost_elements(:,:,3)= cost_elements(:,:,3)./10;
             f=figure;   
%             h=bar3(cost_elements(:,:,2));
%             x=get(h,'Xdata');
%             y=get(h,'Ydata');
           
%             zlabel('infrastructure cost'); 
            ratios=[1 5  20 50 100 200 400 800 1600 3200];
            T= 5;
            for T=1:10
            data=cost_elements(:,T,:);
            data= permute( data,[1,3,2]);
            data=data(1:8,:); 
            data(:,2)= data(:,2)./10;
            h=bar(data);
            axis([.5 8.5 0 900]);
              hXLabel=xlabel('trashing cost coefficient index'); 
             hYLabel=ylabel('horizon'); 
             o.fix_font(hXLabel, hYLabel);
            %hLegend=legend('c_{SLA}',  'c_{infr}',  'c_{trsh}','Location','NorthOutside');
            %set(hLegend,'Interpreter','Latex');
            %set([hLegend, gca],'FontSize', 24 );
            x=get(h,'Xdata');
            y=get(h,'Ydata');
            for xc=1:size(data,1)
                %text(x{2}(xc),y{2}(xc),['y=' num2str(ratios(xc))]);
                text(xc,max(data(xc,:))+20,[num2str(ratios(xc))]);
            end
            r_sla=5;
            r_resource=1;
            objective='abs';
    
            set(gca, ...
              'Box'         , 'off'     , ...
              'TickDir'     , 'out'     , ...
              'TickLength'  , [.02 .02] , ...
              'XMinorTick'  , 'on'      , ...
              'YMinorTick'  , 'on'      , ...
              'XGrid'       , 'on'      , ...              
              'YGrid'       , 'on'      , ...
              'XColor'      , [.3 .3 .3], ...
              'YColor'      , [.3 .3 .3], ... %'XTick'       , 0:20:9, ...
              'LineWidth'   , 2         ,...
              'YTick'       , 0:100:1000, ...
              'LineWidth'   , 2         );

            UtilityLib.print_figure(f,9,7,sprintf('figure/trade_off_mpc%d_%d_%d_%s',...
                T,r_sla, r_resource,objective));
            end 
        
        end

        function cost_elements=data_trade_off_and_horizon(o)

            cost_elements(:,1,:)=[...
             23537.72, 25.16,  0.00
             23537.72, 29.17,  0.00
             23537.72, 104.64,  0.00
             23537.72, 92.09,  0.00
             23537.72, 1446.45,  0.00
             23537.72, 1387.76,  0.00
             23534.35, 2801.14,  0.00
            NaN NaN NaN
            NaN NaN NaN
            NaN NaN NaN
            ];

            cost_elements(:,2,:)=[...
              0.95, 5522.21, 768.00
              0.95, 5685.81, 788.00
              0.95, 8301.54, 95.00
              0.95, 8869.48, 84.00
              0.95, 8994.39, 84.00
             59.11, 8693.31, 84.00
             23537.72, 25.16,  0.00
             23537.72, 276.56,  0.00
             23537.72, 528.34,  0.00
             23537.72, 620.69,  0.00
            ];
            cost_elements(:,3,:)=[...  
              0.95, 5548.36, 746.00
              0.95, 5655.86, 744.00
              0.95, 5902.47, 522.00
              0.95, 8295.46, 95.00
              0.95, 8869.26, 84.00
              0.95, 8993.65, 84.00
             414.52, 7696.79, 86.00
             23537.71, 489.91,  0.00
             23537.72, 213.59,  0.00
             23537.72, 377.39,  0.00
            ];
             cost_elements(:,4,:)=[...
              0.95, 5475.35, 732.00
              0.95, 5638.82, 726.00
              0.95, 5895.17, 456.00
              0.95, 8287.38, 95.00
              0.95, 8868.27, 84.00
              0.95, 8992.85, 84.00
             80.80, 8530.16, 84.00
             788.18, 7151.66, 72.00
             23537.72, 134.10,  0.00
             23537.72, 154.96,  0.00
            ];
            cost_elements(:,5,:)=[...  
              0.95, 5501.41, 744.00
              0.95, 5548.81, 740.00
              0.95, 5812.95, 392.00
              0.95, 8281.31, 95.00
              0.95, 8283.71, 95.00
              0.95, 8868.81, 84.00
             37.18, 8739.79, 84.00
             708.75, 7240.73, 70.00
             23537.44, 3421.67,  0.00
             23537.72, 92.09,  0.00
            ];
            cost_elements(:,6,:)=[...
              0.95, 5505.86, 778.00
              0.95, 5501.60, 778.00
              0.95, 5921.48, 384.00
              0.95, 6122.38, 333.00
              0.95, 8275.38, 95.00
              0.95, 8867.83, 84.00
             64.99, 8886.67, 84.00
             182.68, 8174.77, 79.00
             1020.06, 7000.29, 70.00
             23537.70, 880.57,  0.00
             ];
            cost_elements(:,7,:)=[...
              0.95, 5441.28, 762.00
              0.95, 5510.17, 764.00
              0.95, 5914.93, 410.00
              0.95, 6054.90, 321.00
              0.95, 8269.29, 95.00
              0.95, 8805.65, 84.00
             64.20, 8877.96, 84.00
             181.89, 8058.62, 79.00
             1006.38, 7079.31, 70.00
             23537.66, 1434.08,  0.00
            ];
            cost_elements(:,8,:)=[...
              0.95, 5467.16, 764.00
              0.95, 5475.77, 764.00
              0.95, 5919.18, 382.00
              0.95, 6211.41, 305.00
              0.95, 8263.35, 95.00
              0.95, 8804.19, 84.00
             54.73, 8880.02, 84.00
             169.74, 8140.96, 79.00
             971.60, 6858.41, 70.00
             23537.30, 3454.86,  0.00
            ];
            cost_elements(:,9,:)=[... 
              0.95, 5527.42, 766.00
              0.95, 5480.00, 768.00
              0.95, 5880.65, 372.00
              0.95, 6142.84, 298.00
              0.95, 8257.55, 95.00
              9.46, 8137.85, 98.00
             53.94, 8623.44, 84.00
             155.31, 8218.22, 79.00
             840.95, 7167.06, 70.00
             23537.49, 3471.98,  0.00
            ];
            cost_elements(:,10,:)=[...
              0.95, 5458.55, 774.00
              0.95, 5544.66, 774.00
              0.95, 5850.74, 382.00
              0.95, 6204.77, 288.00
              0.95, 8251.47, 95.00
              9.46, 8131.76, 98.00
             53.94, 8622.48, 84.00
             152.70, 8233.91, 79.00
             812.72, 7152.88, 70.00
             23537.55, 3346.18,  0.00
            ];
        end

    end
    
end


classdef workload
    properties
                header = {'month' 'day' 'hour' 'minute' 'reqnum'};
                m=struct('month', 1 ,'day', 2,'hour', 3,'minute', 4,'reqnum', 5)
    end
    
    methods
        function res=cum_and(obj,m)
            [i j]=size(m);
            res=ones(i,1);
            for k =1:j
                res = res & m(:,k); 
            end                
        end
        
        function res=cum_or(obj,m)
            [i j]=size(m);
            res=zeros(i,1);
            for k =1:j
                res = res | m(:,k); 
            end                
        end
        
        function num=count(obj,m)
            num=size(m,1);
        end
        
        function res=select(obj,m,col_gr,val)
            res = [];
            for k=1:size(m,1)
                if obj.cum_and(m(k,col_gr)==val)
                    res = [res ; m(k,:)];
                end
            end
        end
        
        function [table]=project(obj,m,col_gr)
            B = sortrows(m,col_gr);
            sortedbuckets=B(:,col_gr);
            table = unique(sortedbuckets,'rows');
        end
        
        function [table]=groupby(obj,m,col_gr,col_val)
            B = sortrows(m,col_gr);
            sortedbuckets=B(:,col_gr);
            n=size(sortedbuckets,1);
            distinctbuckets(:,1:n) = 1;
            distinctbuckets(2:n) = obj.cum_or( sortedbuckets(1:n-1,:) ~= sortedbuckets(2:n,:) );
            uniquebuckets=cumsum(distinctbuckets);
            total=accumarray(uniquebuckets',B(:,col_val));
            table = [unique(sortedbuckets,'rows') total];
        end
        
        function [data]=load_data(obj,file)     
            fid = fopen(file);
            tline = fgetl(fid);
            data = [];
            while (ischar(tline))
                [tok mat] = regexp(tline,'(\d+)?-(\d+)?-(\d+)?-(\d+)?\t(\d+)?','tokens', 'match');
                if (~(isempty(tok)))       
                    data = [ data ; str2num(char(tok{:}))'];
                end    
                tline = fgetl(fid);
            end
            fclose(fid);
            data = sortrows(data, [1 2 3 4]);
        end

        function ret=get_workload(obj,file,begining_min,end_min)
                m= obj.load_data(file);                
                header = {'month' 'day' 'hour' 'minute' 'reqnum'};
                lab=struct('month', 1 ,'day', 2,'hour', 3,'minute', 4,'reqnum', 5);
                between = (begining_min : end_min) ;
                ret=m(between ,lab.reqnum);                
        end
        
        function plot_sample(obj)
            workload=obj.get_workload('data/day42_per_min.txt', 1, 23 *60);
            u=UtilityLib();
            handle = figure;
            plot(workload);
            title('Number of users over time');
            xlabel('Time(minutes)');
            ylabel('Number of users');
        end
        
        function plotit(obj)           
                hours = 24;
                m= 1.40*obj.load_data('data/day42_per_min.txt');
                
                header = {'month' 'day' 'hour' 'minute' 'reqnum'};
                lab=struct('month', 1 ,'day', 2,'hour', 3,'minute', 4,'reqnum', 5)
                u=UtilityLib();
                handle = figure;
                
                subplot(3,hours,1:2*hours) 
                between = (1 : 23 *60) ;
                plot(m(between ,lab.reqnum));                
                title('Number of users over time');
                xlabel('Time(minutes)');
                ylabel('Number of users');
                %axis([0 length(res) 0 6])

                pl_num = 2*hours+1;
                coefs = [];
                sols = [];                                
                temp = sortrows(obj.project(m,(1:3)));
                for k=   temp((1:hours),:)'
                    subplot(3,hours,pl_num); 
                    pl_num = pl_num +1;                    
                     res = sortrows(obj.select(m,(1:3),k'),(1:4));
                     %res = sortrows(res,(1:4));   
                     coef = corrcoef(res(:,4),res(:,5));
                     coefs = [coefs coef(1,2)];                     
                     
                     X1 = [ones(size(res,1),1) res(:,4)];
                     y = res(:,5);
                     sol = (X1\y);
                     sols = [sols sol(2,1)]; %(1,1) is the 1s coef
                     line([0;1],[sol(1,1); sol(1,1)+sol(2,1)],[1;1]);
                     axis([0 1 sol(1,1)-1 sol(1,1)+1])
                      set(gca,'YTick',[]); set(gca,'YColor','w');
                      set(gca,'XTick',[]); set(gca,'XColor','w');
                end                

                %set(gca,'Visible','off')
                %u.print_figure(handle,9,7,strcat('./result/figure/',filestr,'/array_instance_number_time'));
        end %plotit 

        % draw_result().get_raw_form('day41-raw.txt','day41_per_minute.txt' )
        function get_raw_form(obj,infile,outfile)           
                hours = 24;
                m= obj.load_data(infile);
                header = {'month' 'day' 'hour' 'minute' 'reqnum'};
                lab=struct('month', 1 ,'day', 2,'hour', 3,'minute', 4,'reqnum', 5)
                
                fid = fopen(outfile,'w');
                between = (3 * 60 : 15 *60) ;              
               fprintf(fid, '%14.0d', round(1.40*m(between,lab.reqnum))); 
               fclose(fid)
        end
   
             % column vector X and Y
        function sol = regress(obj,X,Y)                                                   
                     X1 = [ones(size(X,1),1) X];
                     Y2=Y;
                     windowSize = 5;
                     Y2(isnan(Y))=0;
                     Y2 = filter(ones(1,windowSize)/windowSize,1,Y2')';
                     Y(isnan(Y))=Y2(isnan(Y));    
                     sol = (X1\Y);
                     % sols = [sols sol(2,1)]; %(1,1) is the 1s coef                    
        end
        
      function ret_=draw_rbe_metric_with_regression(obj,a_result)
           u=UtilityLib();    
           delta = 15; %each 15 minutes            
           color_num =1;

            % a_result = a_result(1:219); % limit it to 220 samples
            a_result(a_result<0) = NaN;
             plot(a_result, u.color{color_num});
             tans_ = [];
             sample_per_reg= 15;
             for x1=1:delta:length(a_result)-delta
                x2=x1+delta-1;
                d=obj.regress((1:delta)',a_result(x1:x2));
                intercept = d(1,1);
                tan = d(2,1);                                    
                line([x1; x2],[tan*1+intercept; tan*delta+intercept],...     
                        'Color',u.color{color_num},...
                        'LineStyle',u.lnstyle{color_num},...
                        'LineWidth',2,...
                        'Marker','s',...
                        'MarkerEdgeColor','k',...
                        'MarkerFaceColor','g',...
                        'MarkerSize',6);         
                    tans_=[tans_ tan];
             end
             hold on;
            color_num = color_num+1;
            atan(tans_)*180/pi;
            axis([0 size(a_result,1)  250 max(a_result)]) %hack
            ret_=a_result';
            
        
            title('#Users');
            xlabel('Time');
            ylabel('Number of Users');           
          %  u.print_figure(handle,9,7,strcat(config().log_files.result_dir,'rbe-regression-',metric));
      end

    function ret_=draw_rbe_metric(obj,a_result)
           u=UtilityLib();  
           delta = 15; %each 15 minutes            
           color_num =5;

            % a_result = a_result(1:219); % limit it to 220 samples
            a_result(a_result<0) = NaN;
             plot( (1:size(a_result,1))*250/size(a_result,1),...
                    a_result, 'm');
             tans_ = [];
             sample_per_reg= 15;
             hold on;
            color_num = color_num+1;
            atan(tans_)*180/pi;
          %  axis([0 size(a_result,1)  250 max(a_result)]) %hack
            ret_=a_result';
            
        
            %  title('#Users');
            xlabel('Time');
            ylabel('Clients');           
          %  u.print_figure(handle,9,7,strcat(config().log_files.result_dir,'rbe-regression-',metric));
      end

    
        function handle=draw_raw_form(obj,file)
                fid = fopen(file);
                tline = fgetl(fid);
                data = [];
                data = str2num(tline)';
                fclose(fid);
                
                u=UtilityLib();
                handle = figure;
                
%                 plot(data);                
%                 title('Number of users over time');
%                 xlabel('Time(minutes)');
%                 ylabel('Number of users');
                
                obj.draw_rbe_metric(data(1:720));

%                 pl_num = 2*hours+1;
%                 coefs = [];
%                 sols = [];                                
%                 temp = sortrows(obj.project(m,(1:3)));
%                 for k=   temp((1:hours),:)'
%                     subplot(3,hours,pl_num); 
%                     pl_num = pl_num +1;                    
%                      res = sortrows(obj.select(m,(1:3),k'),(1:4));
%                      %res = sortrows(res,(1:4));   
%                      coef = corrcoef(res(:,4),res(:,5));
%                      coefs = [coefs coef(1,2)];                     
%                      
%                      X1 = [ones(size(res,1),1) res(:,4)];
%                      y = res(:,5);
%                      sol = (X1\y);
%                      sols = [sols sol(2,1)]; %(1,1) is the 1s coef
%                      line([0;1],[sol(1,1); sol(1,1)+sol(2,1)],[1;1]);
%                      axis([0 1 sol(1,1)-1 sol(1,1)+1])
%                      set(gca,'YTick',[]); set(gca,'YColor','w');
%                      set(gca,'XTick',[]); set(gca,'XColor','w');            
%                 end
        end
            
        function test_plot(obj)
            u=UtilityLib();
            handle=obj.draw_raw_form( 'data/day21-raw.txt')
%            u.print_figure(handle,9,7,strcat(config().fig4log_files.result_dir,'img1')); 
%             handle=obj.draw_raw_form('../../bradCNSMData/fig6/day43-raw-partial-mult1.40.txt');
%             u.print_figure(handle,9,7,strcat(config().fig6log_files.result_dir,'img1')); 
        end
        
    end %static methods
end %class



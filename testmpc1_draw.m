classdef testmpc1_draw
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
            % ----------------------------------------------------------------------------
        %  Methods that are used to show the state of placement
        % ----------------------------------------------------------------------------
        
        function ret=print_lat_mat(o, matrix)
            digits(4);
            s = sym(matrix,'d');
            v = vpa(s,5);
            ret=latex(v)
        end
        
        function h=draw_thetaV(o,thetaV,win,T)
            h=figure;
            [x,y,z] = meshgrid( 1:o.num_host ,1:win,1:o.num_service );
            thetaVcut=thetaV(:,:,1:T);
            thetaVcutDiff=thetaV(:,:,2:T+1)-thetaV(:,:,1:T);
            thetaVcut= bsxfun(@rdivide, thetaVcut, max(thetaVcut,[],3));
            thetaVcut(isnan(thetaVcut))=0;
            thetaVcut(thetaVcutDiff==0)=0;
            v= permute(thetaVcut(:,:,[1:win]),[3,1,2]);
            xslice = []; yslice =horzcat([1:2:win]);  zslice =  []; % linspace(1,10,100); linspace(10,T,10)
            slice(x,y,z,v,xslice,yslice,zslice);
            alpha(0.5);
            colormap(gray)
            colormap(flipud(colormap))
            daspect([1,1,1])
            camva(5); camproj perspective;
            campos([-200,-20,56])
            set(findobj(gca,'Type','Surface'),'EdgeColor','none');
            set(gcf,'Color',[.5,.5,.5],'Renderer','zbuffer');
        end
        
        function fig=draw_thetaV2(o,thetaV,T,trsh)
            fig=figure;
            [h s]=find(abs((thetaV(:,:,T)-thetaV(:,:,1))./thetaV(:,:,T))>trsh);
            L={};
            ColorSet = varycolor(size(h,1) );
            set(gca, 'ColorOrder', ColorSet);
            hold all;
            for i=1:size(h,1)
                plot(1:T, permute(thetaV( h(i) , s(i),1:T),[3,1,2]) ,  'Color', ColorSet(i,:) );
                L{i}=sprintf('host%d,serv%d',h(i),s(i));
            end
             hleg1 = legend(L(:));
            title('Transition of Service Placements');
            xlabel('Time'); ylabel('Portion of host multiplicity');
        end
        
        function fig=draw_hosts(o,U_nfm,T)
            fig=figure;
            st = {'-r', '-g',  '--b',  '-.c', '-.m',':k'};
            H={};
            for h=1:o.num_host
                plot(1:T,U_nfm(h,:) ,st{h});  hold on;
                H{h} = sprintf('host%d',h);
            end
            %  hleg1 = legend('host1','host2' ,'host3','host4','host5','host6');
            legend(H(:));
            title('Host Utilizations');
            xlabel('Time'); ylabel('Utilization');
        end
        
        function h=draw_services(o,thetaV,T)
            h=figure;
            plot(permute(sum(thetaV(:,:,1:100),1),[3,2,1]));
            S={};
            for s=1:o.num_service
                S{s} = sprintf('service%d',s);
            end
            %hleg1 = legend('service1','service2' ,'service3','service4','service5','service6',...
            %   'service7','service9' ,'service10','service11','service12','service13','service14' );
            legend(S(:))
            title('Service Resource Association');
            xlabel('Time'); ylabel('Service Resource Association');
        end
        
        function h=draw_classes_throughput(o,X_nfm,X_target,T)
            h=figure;
            plot(1:T,X_nfm(1,:),'-b', 1:T,X_target(1,:),'-.b',  1:T, X_nfm(2,:) ,'-r', 1:T,X_target(2,:),'-.r' );
            hleg1 = legend('X_{c1}', 'X_{c1}^{SLA}', 'X_{c2}', 'X_{c2}^{SLA}');
            title('Class Througput');
            xlabel('Time'); ylabel('Througput');
        end
        
        function h=draw_classes_rt(o,X_nfm,X_target,N,N_t,Z,T)
            h=figure;
            r = bsxfun(@minus, N_t./X_nfm, Z);
            r_target = bsxfun(@minus, N_t./X_target, Z);    %  X_target./N - Z;
            C={};
            ColorSet = varycolor( o.num_classes );
            hold all;
            for c=1:o.num_classes
                plot(1:T,r,'-');
                plot(1:T,r_target,'--');
                C{c} = sprintf('R_{c%d}',c);
            end
            legend(C(:));
            title('Class Response Time');
            xlabel('Time'); ylabel('Response Time');
        end
        
        
        function draw_task_deployment(o, varargin)
            filename = 'test_graph1.dot';
            
            fid = fopen(filename, 'w');
            fprintf(fid, 'digraph G {\n');
            
            for t=1:o.num_container
                [~,s]=find(o.beta_deployed(t,:)>0);
                fprintf(fid, 'T%d [shape=record, label="{{',t);
                for s_i=1:size(s,2);
                    fprintf(fid,'<s%d>s%d', s(s_i),s(s_i));
                    if s_i ~=  size(s,2)
                        fprintf(fid,'|');
                    end
                end
                fprintf(fid, '}}"];\n');
                fprintf(fid, 'T%d -> T%d [label = "T%d" color = transparent]; \n',t,t,t );
            end
            
            for i=1:o.num_host
                fprintf(fid, 'H%d;\n', i  );
                for j=1:o.num_container
                    if o.pl_delta(i,j,1)>0
                        fprintf(fid, 'T%d -> H%d; \n',j,i );
                    end
                end
            end
            
            fprintf(fid, '}\n');
        end
        
        % o= assignment_test_flow_surr()
        % o.draw_service_callgraph()
        % varargin{1} is set to true whenever you want show the
        % deployment on the hosts in a specific timestamp.  the
        % timestamps stored in varargin{2}
        function draw_service_callgraph(o, varargin)
            %   'filename'  -  if omitted, writes to 'tmp.dot'
            o.setup_hosts();
            o.setup_services();
            N=[250 100 200
                100   250 50];
            Z=[1 1]';
            RT_sla=[0.146  0.267]';
            
            i=1;
            filename = 'test_graph1.dot';
            arctxt = '->';
            labeltxt = '[label="%s"]';
            
            fid = fopen(filename, 'w');
            fprintf(fid, 'digraph G {\n');
            %     fprintf(fid, 'node [style=rounded];\n splines=false;');
            fprintf(fid, ' splines=false;');
            
            
            for c=1:o.num_classes
                % fprintf(fid, 'C%d [shape=record, label="{{<c1>Z=..., N=...}| C%d}"];\n',c,c);
                fprintf(fid, 'C%d [shape=record, label="{{<c1>Z=%d, N=%s}}"];\n',c ,Z(c),vect2str(N(c,:), 'formatstring', '%d') );
            end
            
            for t=1:o.num_container
                [~,s]=find(o.beta_deployed(t,:)>0);
                fprintf(fid, 'T%d [shape=record, label="{{',t);
                for s_i=1:size(s,2);
                    fprintf(fid,'<s%d>s%d', s(s_i),s(s_i));
                    if s_i ~=  size(s,2)
                        fprintf(fid,'|');
                    end
                end
                fprintf(fid, '}}"];\n');
                fprintf(fid, 'T%d -> T%d [label = "T%d" color = transparent]; \n',t,t,t );
                
                %  fprintf(fid, '}| T%d}"];\n',t);
            end
            
            for i=1:o.num_classes
                for j=1:o.num_service
                    inv = o.call_graph(i,j+o.num_classes);
                    if  inv>0
                        if varargin{1}==true
                            fprintf(fid,'C%d:c1 ->T%d:s%d; \n',i, o.get_container_of(j), j);
                        else
                            fprintf(fid,'C%d:c1 ->T%d:s%d[label="(%.1f)"]; \n',i, o.get_container_of(j), j, inv);
                        end
                    end
                end
            end
            
            for i=1:o.num_service
                for j=1:o.num_service
                    inv = o.call_graph(i+o.num_classes,j+o.num_classes);
                    if  inv>0
                        if varargin{1}==true
                            fprintf(fid,'T%d:s%d ->T%d:s%d; \n ', o.get_container_of(i), i, o.get_container_of(j), j );
                        else
                            fprintf(fid,'T%d:s%d ->T%d:s%d [label="(%.1f)"]; \n ', o.get_container_of(i), i, o.get_container_of(j), j , inv);
                        end
                    end
                end
            end
            
            timestamp=varargin{2};
            if varargin{1}==true
                for i=1:o.num_host
                    fprintf(fid, 'H%d;\n', i  );
                    for j=1:o.num_container
                        % this is the case when the container is just
                        % added to host and should be drawn in red
                        % link
                        %
                        if (timestamp>1 && o.pl_delta(i,j,timestamp-1)>0)
                            fprintf(fid, 'T%d -> H%d [style=dashed, dir=none,color=blue,penwidth=4]; \n',j ,i );
                        elseif  (timestamp>1 && o.pl_delta(i,j,timestamp-1)<0)
                            fprintf(fid, 'T%d -> H%d [style=dashed, dir=none,color=red,penwidth=4]; \n',j ,i );
                        elseif  (o.placement(i,j, timestamp )>0)
                            fprintf(fid, 'T%d -> H%d [style=dashed, dir=none]; \n',j,i );
                        else
                            %do nothing
                        end
                    end
                end
            end %if
            
            %                      if timestamp>1
            %                          [i j]=find(o.pl_delta(:,:,timestamp-1)~=0);
            %                          c=o.pl_delta(sub2ind(size(o.pl_delta),  i, j, ones(size(i,1),1)*2 ));
            %                          for k=1:size(i,1)
            %                              if (o.pl_delta(i(k),j(k),timestamp-1)==1 )
            %                                  % 'added to';
            %                                  fprintf(fid, 'T%d -> H%d [style=dashed, dir=none,color=red,penwidth=4]; \n',j(k) ,i(k) );
            %                              else
            %                                  % 'removed from';
            %                                  % fprintf(fid, 'T%d -> H%d [style=dashed, dir=none,color=blue,penwidth=4]; \n',j(k) ,i(k) );
            %                              end
            %                          end
            %                      end
            %
            
            
            
            fprintf(fid, '}\n');
            
            
        end
        
             % o= assignment_test_flow_surr()
        % o.draw_service_callgraph()
        % varargin{1} is set to true whenever you want show the
        % deployment on the hosts in a specific timestamp.  the
        % timestamps stored in varargin{2}
        function draw_service_callgraph_nocontainer(o, N, Z, RT_sla, is_deployment, t_)
    %        o.setup_hosts();
    %        o.setup_services();
%             N=[250 100 200
%                 100   250 50];
%             Z=[1 1]';
%             RT_sla=[0.146  0.267]';
    
            %   'filename'  -  if omitted, writes to 'tmp.dot'  
            i=1;
            filename = 'test_graph1.dot';
            arctxt = '->';
            labeltxt = '[label="%s"]';
            
            fid = fopen(filename, 'w');
            fprintf(fid, 'digraph G {\n');
            %     fprintf(fid, 'node [style=rounded];\n splines=false;');
            fprintf(fid, ' splines=false;\n  node [shape="box", penwidth = 2];');
            
            % draws classes
            for c=1:o.num_classes
                % fprintf(fid, 'C%d [shape=record, label="{{<c1>Z=..., N=...}| C%d}"];\n',c,c);
                fprintf(fid, 'C%d [shape=record, label="{{<c1>Z=%d, N=%s}}"];\n',c ,Z(c),vect2str(N(c,:), 'formatstring', '%d') );
            end
            
            % services 
            for s=1:o.num_service 
                fprintf(fid, 's%d [style="rounded"] \n',s );
            end
            
            % arrows from classes to services 
            for i=1:o.num_classes
                for j=1:o.num_service
                    inv = o.call_graph(i,j+o.num_classes);
                    if  inv>0
                         if is_deployment==true
                            fprintf(fid,'C%d:c1 ->s%d; \n',i, j);
                         else
                            fprintf(fid,'C%d:c1 -> s%d[label="(%.1f)"]; \n',i,  j, inv);
                        end
                    end
                end
            end
          
            % arrows from services to services 
            for i=1:o.num_service
                for j=1:o.num_service
                    inv = o.call_graph(i+o.num_classes,j+o.num_classes);
                    if  inv>0
                          if is_deployment==true
                            fprintf(fid,'s%d ->s%d; \n ', i,  j );
                        else
                            fprintf(fid,'s%d ->s%d [label="(%.1f)"]; \n ', i, j , inv);
                        end
                    end
                end
            end
            
            % placement arrows from services to hosts  
            if is_deployment==true
                timestamp=t_;
                for i=1:o.num_host
                    fprintf(fid, 'H%d;\n', i  );
                    for j=1:o.num_service
                        % this is the case when the container is just
                        % added to host and should be drawn in red
                        % link
                        %
                        if (timestamp>1 && o.pl_delta(i,j,timestamp-1)>0)
                            fprintf(fid, 's%d -> H%d [style=dashed,  dir=none,color=blue,penwidth=4]; \n',j ,i );
                        elseif  (timestamp>1 && o.pl_delta(i,j,timestamp-1)<0)
                            fprintf(fid, 's%d -> H%d [style=dashed, dir=none,color=red,penwidth=4]; \n',j ,i );
                        elseif  (o.placement(i,j, timestamp )>0)
                            fprintf(fid, 's%d -> H%d [style=dashed, dir=none]; \n',j,i );
                        else
                            %do nothing
                        end
                    end
                end
            end %if      
            fprintf(fid, '}\n');    
        end
   
        
        function print_delta_sequece(o,pl_delta)
            s=size(pl_delta);
            fprintf('\\begin{eqnarray*}\n\\begin{split}\n');
            for t=2:s(3)
                fprintf('\\text{t=%d:} \\\\\n',t);
                [i j]=find(pl_delta(:,:,t)~=0);
                c=pl_delta(sub2ind(size(pl_delta),  i, j, ones(size(i,1),1)*2 ));
                for k=1:size(i,1)
                    if (pl_delta(i(k),j(k),t)==1 )
                        text='added to';
                    else
                        text='removed from';
                    end
                    fprintf('& %d   &  \\text{%s }  &  %d\\\\ \n',j(k) ,text,i(k))
                end
            end
            fprintf('\\end{split}\n\\end{eqnarray*}\n');
            
            %             \begin{eqnarray*}
            %             \begin{split}
            %
            %             \text{t=2:} \\
            %             & 8   &  \text{added to }  & 2\\
            %             & 8  &   \text{added to } & 4 \\
            %
            %             \text{t=3:}  \\
            %              &  3   &  \text{remove from } & 1 \\
            %             &   7  &    \text{added to } & 4 \\
            %              &   8  &    \text{remove from } & 4 \\
            %             \end{split}
            %             \end{eqnarray*}
        end
        
        function print_mat( obj, mat )
            digits(4);
            s = sym(mat,'d');
            v = vpa(s,5);
            latex(v)
        end
        
      
    end
    
end


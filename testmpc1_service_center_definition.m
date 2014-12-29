classdef testmpc1_service_center_definition
    %UNTITLED7 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
            
        % ----------------------------------------------------------------------------
        %  Methods  for  setting of different types of  service center
        % ----------------------------------------------------------------------------
          
       function setup_hosts(o)
            o.num_host=6;
            o.cap=[8 9 9 8 9 6]'; %[6 9 12]';
            o.cap=o.cap+1;
            o.speed_ratio=[1 1.2 0.9 1.1 0.8 1.2]';
            RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',sum(10)));
            % rng(sum(100*clock),'v4')
            o.c=rand(o.num_host,o.num_service)+10;
            o.c;
            o.cost=round(rand(o.num_host,1)*10)+1;
            % c=round(rand(num_host,1)*10);
        end
        
        function setup_services(o)
            o.num_service=14;
            o.num_classes = 2;
            
            call=zeros(o.num_service+o.num_classes,o.num_service+o.num_classes);
            call(1,3)=1; call(1,4)=2;  call(1,5)=1;
            call(2,5)=2; call(2,6)=1;
            call(3,7)=0.8;
            call(4,7)=1;
            call(5,8)=2.5; call(5,9)=1;
            call(6,8)=2; call(6,10)=1; call(6,15)=0.5;
            call(7,11)=1;
            call(8,11)=1; call(8,12)=2;
            call(9,13)=0.5; call(9,14)=1;
            call(11,16)=1.5;
            o.call_graph = call;
            Y=calculate_contact_level(o,call);
            
            
            D=...
                [0.005, 0.004, 0.002, 0.008, ...
                0.001, 0.002, 0.005, 0.008, ...
                0.006, 0.008, 0.004, 0.005, 0.006,...
                0.005];
            
            o.d= Y .* (ones(o.num_classes,1) * D);
            
            
        end % setup_services
        

    end
    
end


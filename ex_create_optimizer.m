


            
function session=ex_create_optimizer(num_host,cap,speed_ratio,num_classes,...
    num_service,r_sla, r_resource, r_trsh,  d)
            o=testmpc1(); 
            o.num_host=num_host; %6
            o.cap=cap'; %[8 9 9 8 9 6]'; 
            o.speed_ratio=speed_ratio'; %[1 1.2 0.9 1.1 0.8 1.2]';
            o.num_service=num_service; %14;
            o.num_classes =num_classes; % 2;
            
            RandStream.setDefaultStream ...
             (RandStream('mt19937ar','seed',10));
            c_orig=rand(num_host,1)+10; 
            c=c_orig*ones(1, num_service);
           o.c=c;
           
      
        
            o.d= d;

            %session=ex_optimizer_session();
            %session.optimizer=o;
            %-- session.thetaV_0 = zeros(o.num_host , o.num_service);
            %session.r=[r_sla r_resource r_trsh];           
            session=struct('optimizer',o,'r',[r_sla r_resource r_trsh]);
end



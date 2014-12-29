        
num_host=3;
num_service=4;
 num_classes = 2; 
d= rand(num_service,num_classes)
theta= rand(num_service,num_host); 
d_csh= permute(  repmat(d, [ 1,1 , num_host] ) , [2,1,3])  .* ... 
 permute( repmat(theta, [1 , 1, num_classes]), [3,1,2]  ); 

 P = zeros(2,3,2,3);
 P(1,1,1,3) = 1;
 P(1,3,1,1) = 1;
 P(2,2,2,3) = 1;
 P(2,3,2,2) = 1;
 V = qncmvisits(P,[3 3]) % reference station is station 3        
 
 
 
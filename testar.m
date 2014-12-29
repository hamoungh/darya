% http://stackoverflow.com/questions/12925443/generating-appropriately-scaled-ar-process-using-matlab-filter

A=[1 -2.7607 3.8106 -2.6535 0.9238]; 
%A = [1 -0.7 -0.25];
% AR(4) coefficients
y=filter(1,A,sqrt(0.2)*randn(100,1)); 
% Filter a white noise input to create AR(4) process
[ar_coeffs,nv] =arburg(y,4)
%compare the results in ar_coeffs to the vector A.
plot(y) 


x=randn(1000,1);
rho=0.50;

% Old while-loop      
y2=x;
I=2;
while I<=length(x);
   y2(I)=10+y2(I-1)*rho+x(I);
   I=I+1;
end;

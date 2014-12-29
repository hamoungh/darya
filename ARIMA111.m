% http://www.mathworks.com/support/solutions/en/data/1-15WF0/?product=ML 
% yin=10+sin([1:100]/10);
% yy=ARIMA111(yin);
% 
function Y = ARIMA111()
A=[1 -2.7607 3.8106 -2.6535 0.9238]; 
%A = [1 -0.7 -0.25];
% AR(4) coefficients
y=filter(1,A,sqrt(0.2)*randn(100,1)); 
% Filter a white noise input to create AR(4) process
[ar_coeffs,nv] =arburg(y,4)
%compare the results in ar_coeffs to the vector A.
plot(y) 

predictions=1;
p=4;
q=0;
%yy = y(1);
y = diff(y);
Spec = garchset('R',p,'M',q,'VarianceModel','Constant');
[EstSpec,EstSE] = garchfit(Spec,y); 
EstSpec
Y = garchsim(EstSpec,50);
%[sigmaForecast,meanForecast] = garchpred(EstSpec,y,predictions);
%pred = cumsum([yy; y; meanForecast])
%pred = pred(end-predictions+1:end);
end

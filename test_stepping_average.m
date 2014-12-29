t = sort((100:999)' + 3*rand(900,1));     % non-uniform time
x = 5*rand(900,1) + 10;             % x(i) is the value at time t(i)

tt = ( floor(t(1)):1*60:ceil(t(end)) )';
int = 1 + floor((t - t(1))/60);
mu = accumarray(int,x,[],@mean)
sd = accumarray(int,x,[],@std);
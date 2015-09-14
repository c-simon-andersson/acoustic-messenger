function y = pulsetr(fun,a, tau,n,width,data)

%fun: function describing the pulse
%a: function parameter
%n: number of samples per pulse
%width: total width of pulse support in symbol intervals, must be multiple
%of 1/n
%data: 

%Change symbol time here
%tau = 1;

%Compute working variables
int = tau/n; len = length(data); tmp = tau*width/2;
t = -tmp:int:tmp; lt = length(t);
%t is basic pulse points
num = n*(width+len-1) + 1;
%Total output width
%Superpose pulses one by one, left to right
y = zeros(1,num); x = feval(fun,a, tau, t);
for k = 1:len,
    tmp = n*(k - 1);
    y = y + [zeros(1,tmp) x*data(k) zeros(1,num-lt-tmp)];
end


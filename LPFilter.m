function y = LPFilter(x,fc)
k = .73;
alpha = (1-k*cos(2*pi*fc)-sqrt(2*k*(1-cos(2*pi*fc))-k^2*sin(2*pi*fc)^2))/(1-k);
y = zeros(size(x));
for i = 1:size(x,1)
    y(i,:) = filtfilt(1-alpha,[1 -alpha],x(i,:));
end
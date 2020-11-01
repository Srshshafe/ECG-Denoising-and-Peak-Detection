function [phase phasepos] = PhaseCalculation(peaks)
phasepos = zeros(1,length(peaks));
I = find(peaks);
for i = 1:length(I)-1;
    m = I(i+1) - I(i);
    phasepos(I(i)+1:I(i+1)) = 2*pi/m : 2*pi/m : 2*pi;
end
m = I(2) - I(1);
L = length(phasepos(1:I(1)));
phasepos(1:I(1)) = 2*pi-(L-1)*2*pi/m:2*pi/m:2*pi;
m = I(end) - I(end-1);
L = length(phasepos(I(end)+1:end));
phasepos(I(end)+1:end) = 2*pi/m:2*pi/m:L*2*pi/m;
phasepos = mod(phasepos,2*pi);
phase = phasepos; 
I = find(phasepos>pi);
phase(I) = phasepos(I) - 2*pi;



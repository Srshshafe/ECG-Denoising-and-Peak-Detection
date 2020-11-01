function [ECGmean,ECGsd,meanPhase] = MeanECGExtraction(x,phase,bins,flag)
meanPhase = zeros(1,bins);
ECGmean = zeros(1,bins);
ECGsd = zeros(1,bins);
I = find( phase>=(pi-pi/bins)  | phase<(-pi+pi/bins) );
if(~isempty(I))
    meanPhase(1) = -pi;
    ECGmean(1) = mean(x(I));
    ECGsd(1) = std(x(I));
else
    meanPhase(1) = 0;
    ECGmean(1) =0;
    ECGsd(1) = -1;
end
for i = 1 : bins-1;
    I = find( phase >= 2*pi*(i-0.5)/bins-pi & phase < 2*pi*(i+0.5)/bins-pi );
    if(~isempty(I))
        meanPhase(i + 1) = mean(phase(I));
        ECGmean(i + 1) = mean(x(I));
        ECGsd(i + 1) = std(x(I));
    else
        meanPhase(i + 1) = 0;
        ECGmean(i + 1) = 0;
        ECGsd(i + 1) = -1;
    end
end
K = find(ECGsd==-1);
for i = 1:length(K),
    switch K(i)
        case 1
            meanPhase(K(i)) = -pi;
            ECGmean(K(i)) = ECGmean(K(i)+1);
            ECGsd(K(i)) = ECGsd(K(i)+1);
        case bins
            meanPhase(K(i)) = pi;
            ECGmean(K(i)) = ECGmean(K(i)-1);
            ECGsd(K(i)) = ECGsd(K(i)-1);
        otherwise
            meanPhase(K(i)) = mean([meanPhase(K(i)-1),meanPhase(K(i)+1)]);
            ECGmean(K(i)) = mean([ECGmean(K(i)-1),ECGmean(K(i)+1)]);
            ECGsd(K(i)) = mean([ECGsd(K(i)-1),ECGsd(K(i)+1)]);
    end
end
if(flag==1)
    ECGmean = ECGmean - mean(ECGmean(1:ceil(length(ECGmean)/10)));
end
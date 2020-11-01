clc
clear all
close all
load('ECG_Sample.mat'); 
data = data(1:15000,2)';
fs = 1000;
t = (0:length(data)-1)/fs;
f = 1.2;
% removing baseline using a low pass filter function
bsline = LPFilter(data,.73/fs);
data1 = data - bsline;
% Creating noise
SNR = 15;
SignalPower = mean(data1.^2);
NoisePower = SignalPower / 10^(SNR/10);
x = data1 + sqrt(NoisePower)*randn(size(data1));
peaks = PeakDetection(x,f/fs);  
I = find(peaks);
[phase phasepos] = PhaseCalculation(peaks);
teta = 0;
pphase = PhaseShifting(phase,teta);  
bins = 250;
[ECGmean,ECGsd,meanphase] = MeanECGExtraction(x,pphase,bins,1);
load ModelParams
%ECGBeatFitter(ECGmean,ECGsd,meanphase,'OptimumParams');
L = (length(OptimumParams)/3);
ai = OptimumParams(1:L); 
bi = OptimumParams(L+1:2*L);
thetai = OptimumParams(2*L+1:3*L); 
[T,X] = ode45(@ECG,(0:0.0001:10),[0.7,0.7,0.0001],[],thetai,ai,bi);
plot3(X(:,1),X(:,2),X(:,3)+sqrt(NoisePower)*randn(size(X(:,3))))
grid on
xlim([-1 1]),ylim([-1 1]) 
N = length(OptimumParams)/3; % number of Gaussian kernels
JJ = find(peaks);
fm = fs./diff(JJ);          % heart-rate frequancy
w = mean(2*pi*fm);          % average heart-rate in rads 
wsd = std(2*pi*fm,1);       % heart-rate standard deviation in rads.
y = [phase ; x];
X0 = [-pi 0]'; 
P0 = [(2*pi)^2 0 ;0 (10*max(abs(x))).^2];
Q = diag( [ (.1*OptimumParams(1:N)).^2 (.05*ones(1,N)).^2 (.05*ones(1,N)).^2 (wsd)^2 , (.05*mean(ECGsd(1:round(length(ECGsd)/10))))^2] );
R = [(w/fs).^2/12 0 ;0 (mean(ECGsd(1:round(length(ECGsd)/10)))).^2];
Wmean = [OptimumParams w 0]';
Vmean = [0 0]';
Inits = [OptimumParams w fs];
InovWlen = ceil(.5*fs);     % innovations monitoring window length
tau = [];                   % Kalman filter forgetting time. tau=[] for no forgetting factor
gamma = 1;                  % observation covariance adaptation-rate. 0<gamma<1 and gamma=1 for no adaptation
RadaptWlen = ceil(fs/2);    % window length for observation covariance adaptation
Xhat = ExtendedKalmanFilter(y,X0,P0,Q,R,Wmean,Vmean,Inits,InovWlen,tau,gamma,RadaptWlen,1);
Xhat = Xhat(2,:);
figure
plot(t,x,'--b'); 
hold on;
plot(t,Xhat','r');
grid;
legend('Noisy','EKF Output');
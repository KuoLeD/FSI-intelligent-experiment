function Coe = DataDeal(data,U,A,f,L,D,M,tinterval,CF_CL,IL_CL)
% This program processes forced-oscillation experimental data
% First, perform filtering on the data
Fs = 1000;  % sampling frequency is 1000 Hz
dt = 1/Fs;
fc = 1;  % set cutoff frequency to 1 Hz
[b, a_CF] = butter(4, fc/(Fs/2), 'low');  % 4th-order low-pass Butterworth filter
data(end,:)=[];
for i=7:10
    % Filter the data
    data(:,i) = filtfilt(b, a_CF, data(:,i));
end
% Data processing -- single cylinder
a4=data(:,4); % dis CF
a5=data(:,5); % dis IL
a6=data(:,6); % vel CF
a7=data(:,7); % F_x1 CF upper end
a8=data(:,8); % F_y1 IL upper end
a9=data(:,9); % F_x2 CF lower end
a10=data(:,10); % F_y2 IL lower end
% Determine the range of the initial segment
N=length(data(:,1)); % total data length
w0=Fs*1/4*tinterval;
w1=Fs*3/4*tinterval;
range0=w0:w1;        % take the latter 1/2 tinterval as the stable-segment range
range0il=8*w0:8*w1;
fx10=a7(range0);
fy10=a8(range0);
fx20=a9(range0il);
fy20=a10(range0il);
dis0_CF=a4(range0);
dis0_IL=a5(range0);

fxm0=mean(fx10+fx20); % mean of the initial segment
fym0=mean(-fy10+fy20);
dism0_CF=mean(dis0_CF);
dism0_IL=mean(dis0_IL);

% Identify the stable-segment range
% Since the data fluctuation is irregular and cannot be identified by peaks/valleys, the stable segment is manually set here
n0=Fs*tinterval; % static sampling for tinterval seconds
n=round((N-n0)/10);
range1=n0+n*2:n0+n*9; % after the test starts, trim the beginning and end; stable segment 20%-80%

fx1mean=mean(a7(range1));
range10 = range1(0.01>=(a7(range1)-fx1mean) & (a7(range1)-fx1mean)>= -0.01);
range1_min = min(range10);
range1_max = max(range10);
fx1=a7(range1_min:range1_max);
fy1=a8(range1_min:range1_max);
fx2=a9(range1_min:range1_max);
fy2=a10(range1_min:range1_max);
dis1_CF=a4(range1_min:range1_max);
dis1_IL=a5(range1_min:range1_max);

fx=fx1+fx2-fxm0;
fy=-fy1+fy2-fym0;
dis_CF=dis1_CF-dism0_CF;
dis_IL=dis1_IL-dism0_IL;

fxm=mean(fx); % mean of the stable segment
fym=mean(fy);
dism_CF=mean(dis_CF);
dism_IL=mean(dis_IL);

% Amplitude after filtering
disp_CF=bpass_CF(dis_CF,dt,0.001,5);
disp_IL=bpass_CF(dis_IL,dt,0.001,5);

rou=1000;
Cdm=fym/(0.5*U*U*L*D*rou);
CLm=fxm/(0.5*U*U*L*D*rou);
if Cdm>3
   YoN=0;
else
   YoN=1;
end
%% fft of sway CF Dis
m=length(disp_CF); % number of samples; fft resolution is Fs/m
ff_sway=Fs*(0:m/2)/m;
ff_sway=ff_sway';
P2=abs(fft((disp_CF))/m);   % normalized fft magnitude corresponding to frequency; two-sided spectrum
P1 = P2(1:round(m/2)+1);
P1(2:end-1)=2*P1(2:end-1);
sway_fft=P1;
[~,p_sway]=max(sway_fft);
f_sway_CF=ff_sway(p_sway);
% figure(3)
% plot(ff_sway(1:floor(m/100)),sway_fft(1:floor(m/100)),'k');
% title('(DISP)','FontName','Times New Roman','FontSize',12)
% ylabel('Mag.(ms)','FontName','Times New Roman','FontSize',12)
% xlabel('Fre. (Hz)','FontName','Times New Roman','FontSize',12)
% set(gca,'FontName','Times New Roman','FontSize',12)
%% FFT of sway IL dis
m=length(disp_IL);
ff_sway=Fs*(0:m/2)/m;
ff_sway=ff_sway';
P2=abs(fft(disp_IL)/m);
P1=P2(1:round(m/2)+1);
P1(2:end-1)=2*P1(2:end-1);
sway_fft=P1;
[~,p_sway]=max(sway_fft);
f_sway_IL=ff_sway(p_sway);

%% phase
f_uper_CF=f_sway_CF*5; % 5 times the dominant frequency of CF displacement
disp_CF=bpass_CF(disp_CF,dt,0,f_uper_CF);
x_CF=disp_CF*0.125/50000; % 0.125/50000 converts pulses to displacement
f_uper_IL=f_sway_IL*5; % 5 times the dominant frequency of IL displacement
disp_IL=bpass_IL(disp_IL,dt,0,f_uper_IL);
y_IL=disp_IL*0.125/50000; % 0.125/50000 converts pulses to displacement
% Assume t is the time vector and x is the cross-flow displacement vector
h = dt; % step size
n = length(x_CF);  % data length
% Cross-flow
a_CF = zeros(n,1); % initialize acceleration vector
% Compute acceleration using the second-order central difference method
for i = 2:n-1
    a_CF(i) = (x_CF(i+1) - 2*x_CF(i) + x_CF(i-1)) / h^2;
end
% Handle boundary points
a_CF(1) = a_CF(2);
a_CF(n) = a_CF(n-1);
% In-line
a_IL = zeros(n,1); % initialize acceleration vector
% Compute acceleration using the second-order central difference method
for i = 2:n-1
    a_IL(i) = (y_IL(i+1) - 2*y_IL(i) + y_IL(i-1)) / h^2;
end
% Handle boundary points
a_IL(1) = a_IL(2);
a_IL(n) = a_IL(n-1);

Fx_ma_CF=M*a_CF; % added-mass inertia force in CF
Fy_ma_IL=M*a_IL; % added-mass inertia force in IL
% Remove mean (i.e., average offset)
Fx=bpass_CF(fx,dt,0.0001,f_uper_CF); % CF
Fy=bpass_CF(fy,dt,0.0001,f_uper_IL); % IL (mean drag removed)
% Subtract inertia force
Fxr=Fx-Fx_ma_CF;
Fyr=Fy-Fy_ma_IL;
% Set vibration frequency band
flow_CF=f_sway_CF*0.9;
fhigh_CF=f_sway_CF*1.1;
flow_IL=f_sway_IL*0.9;
fhigh_IL=f_sway_IL*1.1;
% Obtain displacement and force within the set vibration frequency band (ignore forces at other frequencies, St force)
% disp_CF=bpass_CF(x_CF,dt,flow_CF,fhigh_CF); 
% disp_IL=bpass_CF(y_IL,dt,flow_IL,fhigh_IL);
Fx_main=bpass_CF(Fxr,dt,flow_CF,fhigh_CF);
Fy_main=bpass_CF(Fyr,dt,flow_IL,fhigh_IL);
% Displacement amplitudes in CF and IL
A_CF=rms(x_CF)*sqrt(2)/D;
A_IL=rms(y_IL)*sqrt(2)/D;
% Lift-force amplitude and coefficient amplitude in CF and IL
AFx=rms(Fx_main)*sqrt(2);
CL0_CF=AFx/(0.5*U*U*L*D*rou);
AFy=rms(Fy_main)*sqrt(2);
CL0_IL=AFy/(0.5*U*U*L*D*rou);
% Phase difference between displacement and force in CF and IL directions
disph_CF=hilbert(x_CF);
Fxh=hilbert(Fx_main);
phasex=unwrap(angle(disph_CF));
phaseFx=unwrap(angle(Fxh));
deltaphase_CF=phaseFx-phasex;

disph_IL=hilbert(y_IL);
Fyh=hilbert(Fy_main);
phasey=unwrap(angle(disph_IL));
phaseFy=unwrap(angle(Fyh));
deltaphase_IL=phaseFy-phasey;

% Discard end effects at both ends of the phase-difference time history
deltaphase1=deltaphase_CF(length(deltaphase_CF)*0.1:length(deltaphase_CF)*0.9);
deltaphase2=deltaphase_IL(length(deltaphase_IL)*0.1:length(deltaphase_IL)*0.9);
%% Lift and added-mass coefficients in CF direction
phi1=mean(deltaphase1);
CLV_CF=CL0_CF*sin(phi1);
Cla_CF=-CL0_CF*cos(phi1);
f0non_CF=f_sway_CF*D/U;
Cma_CF=-1/(2*pi^3)*Cla_CF/(A_CF*f0non_CF^2);
% Lift and added-mass coefficients in IL direction
phi2=mean(deltaphase2);
CLV_IL=CL0_IL*sin(phi2);
Cla_IL=-CL0_IL*cos(phi2);
f0non_IL=f_sway_IL*D/U;
Cma_IL=-1/(2*pi^3)*Cla_IL/(A_IL*f0non_IL^2);
% if Cma>4
%    YoN=0;
% else
%    YoN=1;
% end
% if CLV_CF>2
%    YoN=0;
% else
%    YoN=1;
% end
if CF_CL==0
    CLV_CF=0;
    Cma_CF=0;    
end
if IL_CL==0
    CLV_IL=0;
    Cma_IL=0;
end
Coe=[CLV_CF,Cma_CF,CLV_IL,Cma_IL,Cdm];
return

%% Silicon Ellipsometry at azimuthal angle 45 degrees
%% Presetting
clc;
clear all;  %clear the memory
close all;  %close the previous plots and or data matrixes

%% Datas 

c=3E8; %m/s

TETA=30; %degrees - Change accordingly
QHI=-35; %degrees - Change accordingly

%Create strings
PATH1 = strcat('P_QHI=',string(QHI),',','TETA=',string(TETA),'.txt');
PATH2 = strcat('S_QHI=',string(QHI),',','TETA=',string(TETA),'.txt');

%% Importing datas
DATAP = importdata(PATH1); %P 
DATAS= importdata(PATH2); %S 

Time = DATAP(:,1);%s %Setting first column as time
Time = Time-Time(1);%s %Setting firt value to 0 and subsequent
Time = Time.*1E12; %In ps

P = DATAP(:,2);%a.u. %Setting second column as the amplitude 
S = DATAS(:,2);%a.u. %Setting second column as the amplitude

%% Remove Background
TimeStep = Time(end) - Time(end-1);
Time_index = 5; %ps %We are going to do a mean from 0 up to Time_index
End_Index = floor(Time_index./TimeStep);

Background_1 = mean(P(1:End_Index));
Background_2 = mean(S(1:End_Index));

P = P-Background_1;
S = S-Background_2;

%% Zero Padding - Not used
% TimeStep = Time(end) - Time(end-1);
% Cut_Index1 = floor(14./TimeStep); 
% Cut_Index2 = floor(25./TimeStep); 
% Sample_len = length(DATA1P(:,2));
% P = [zeros(Cut_Index1-1,1);DATA1P(Cut_Index1:Cut_Index2,2);zeros(Sample_len-(Cut_Index2),1)];
% S = [zeros(Cut_Index1-1,1);DATA1S(Cut_Index1:Cut_Index2,2);zeros(Sample_len-(Cut_Index2),1)];

% TimeStep = Time(end) - Time(end-1);
% Cut_Index2 = floor(25./TimeStep); 
% Sample_len = length(DATA1P(:,2));
% P = [DATA1P(1:Cut_Index2,2);zeros(Sample_len-(Cut_Index2),1)];
% S = [DATA1S(1:Cut_Index2,2);zeros(Sample_len-(Cut_Index2),1)];

%% Figure 1 - E [a.u.] vs Time [ps]
    figure(1)
    plot(Time,P); %Plot E of the P vs Time (120ps)
    hold on
    plot(Time,S); %Plot E of the S vs Time (120ps)
    hold off
    
    title_string = strcat('THz Time Domain Signal (QHI,TETA)=', '(' ,string(QHI), ',', string(TETA),')');
    title(title_string);
    xlabel('Time [ps]')
    ylabel('E Field [a.u]')
    legend({'P','S'},'Location','southeast')
    
    saveas(gca, title_string);

%% Variable change
    x = Time;
    y = P;  
    z = S;

%% FFT 
                          
    f_sampling = 1E12/(x(end)-x(end-1)); %1/Delta_t = Frequency Resolution [THz]
    L=length(y); %Number of elements in the vector
    f =(f_sampling*(0:(L/2))/L);
    f=f';
%     f = linspace(0,f_sampling/2,L/2+1)'; %Column of frequency axis

%% FFT Amplitudes    

    %FFT P POL
    FFT_P_POL = fft(y,L)/L; % Normalized to L
    %AMPLITUDE OF FFT 
    AFFT_P_POL = [abs(FFT_P_POL(1));2*abs(FFT_P_POL(2:L/2+1))];  %(Spalma la potenza) 

    %FFT SAMPLE                      
    FFT_S_POL = fft(z,L)/L;
    %AMPLITUDE OF FFT sample signal
    AFFT_S_POL = [abs(FFT_S_POL(1));2*abs(FFT_S_POL(2:L/2+1))]; 
 
      
%% Figure 2 - Amplitude Spectra dB vs THz
    figure(2)
    
    semilogy(f./1E12,AFFT_P_POL);
    hold on
    semilogy(f./1E12,AFFT_S_POL);
    hold off
    
    xlim([0.1 2]) 
    title_string = strcat('Amplitude Spectra (QHI,TETA)=', '(' ,string(QHI), ',', string(TETA),')');
    title(title_string);
    xlabel('Frequency[THz]')
    ylabel('Amplitude [dB]')
    legend({'P POL','S POL'},'Location','southeast')
    
    saveas(gca, title_string);

%% FFT Phases
    
    X2=FFT_P_POL;
    phase=atan2(imag(X2),real(X2))*180/pi; %Four-quadrant arctan in degrees
    phase_single=phase(1:L/2+1); %turns the plot into a one-sided plot
    
    %phase_single(2:end-1)=phase_single(2:end-1);
    
    %Unwrapping the phase
    Phase_P_POL= (180/pi)*unwrap(phase_single*pi/180); %PHASE OF FFT in degrees
    
  
    X3=FFT_S_POL;
    phase=atan2(imag(X3),real(X3))*180/pi;
    phase_single=phase(1:L/2+1); %turns the plot into a one-sided plot
    
    %phase_single(2:end-1)=phase_single(2:end-1);

    %Unwrapping the phase
    Phase_S_POL= (180/pi)*unwrap(phase_single*pi/180); %PHASE OF FFT in degrees


%% Figure 3 - Phase [degrees] vs Frequency [THz]
    figure(3)
    
    plot(f./1E12,Phase_P_POL);
    hold on 
    plot(f./1E12,Phase_S_POL); 
    hold off
    
    xlim([0.1 2]) 
    title_string = strcat('FFT Phases (QHI,TETA)=', '(' ,string(QHI), ',', string(TETA),')');
    title(title_string);
    xlabel('Frequency[THz]')
    ylabel('phase')
    legend({'P_POL','S_POL'},'Location','southeast')

   saveas(gca, title_string);

%% EXTRACTING Complex Parameters

w=2.*pi*f; %Vector of the angular frequencies
 
Tan_PHI=AFFT_P_POL./AFFT_S_POL;%Amplitude Ratio
Delta_phase = Phase_P_POL-Phase_S_POL; %!n.b. This is in degree

%We suppose, since we have already calculated it in transmissiont that
%@1THz we have n = 3.4;

n_r = 3.4;

%Analytical calculations
r_p = (sqrt(n_r^2-(sin(TETA*pi/180))^2)-n_r^2*cos(TETA*pi/180))/(sqrt(n_r^2-(sin(TETA*pi/180))^2)+n_r^2*cos(TETA*pi/180));
r_s = (cos(TETA*pi/180)-sqrt(n_r^2-(sin(TETA*pi/180))^2))/(cos(TETA*pi/180)+sqrt(n_r^2-(sin(TETA*pi/180))^2));

Ratio_Analytical = r_p/r_s;

%Ratio from data @1THz
Cut_Index = floor(1E12/(f(2)-f(1)));
Ratio_experimental = Tan_PHI(Cut_Index);

%To Save the values
Values = [QHI,TETA,r_p,r_s,Ratio_Analytical,Ratio_experimental];
dlmwrite('Check Compensation.txt',Values,'-append');
  
%% EXTRACTING DIELECTRIC CONSTANT INDEX OF REFRACTIONS

%Real Part 

SIN2 = (sin(TETA.*pi./180)).^2;
TAN2 = (tan(TETA.*pi./180)).^2;
COS2PHI2 = (cos(2.*atan(1./Tan_PHI))).^2;
SIN2PHI2 = (sin(2.*atan(1./Tan_PHI))).^2;
SIN2PHI = (sin(2.*atan(1./Tan_PHI)));
SIN4PHI = (sin(4.*atan(1./Tan_PHI)));
COSDELTA = cos(Delta_phase.*pi./180);
SINDELTA = sin(Delta_phase.*pi./180);

epsilon_r = SIN2.*(1+(TAN2.*(COS2PHI2-SIN2PHI2.*SINDELTA.^2)))./((1+SIN2PHI.*COSDELTA).^2);
epsilon_r_mean = epsilon_r(1:Cut_Index);
epsilon_r_mean = ones(L/2+1,1).*mean(epsilon_r_mean);

epsilon_i = SIN2.*((TAN2.*SIN4PHI.*SINDELTA)./(1+SIN2PHI.*COSDELTA).^2);
epsilon_i_mean = epsilon_i(1:Cut_Index);
epsilon_i_mean = ones(L/2+1,1).*mean(epsilon_i_mean);

n=sqrt(epsilon_r);
n_mean = n(1:Cut_Index);
n_mean = ones(L/2+1,1).*mean(n_mean);

k=1./(2.*n).*(epsilon_i);
k_mean = k(1:Cut_Index);
k_mean = ones(L/2+1,1).*mean(k_mean);

%% Figure 4 - Index of refractions vs Frequency [THz]
    figure(4)
    
    a = plot(f./1E12,n);
    hold on 
    b = plot(f./1E12,k);
    plot(f./1E12,n_mean);
    plot(f./1E12,k_mean);
    hold off
    
    xlim([0.1 1]) 
    %ylim([-1 15]) 
    
    title_string = strcat('Index of refraction SUBSTRATE (QHI,TETA)=', '(' ,string(QHI), ',', string(TETA),')');
    title(title_string);
    xlabel('Frequency[THz]')
    ylabel('n&k')
    legend([a,b],{'n','k'},'Location','southeast')
    
    saveas(gca, title_string);

%% Figure 5 - Dielectric constant vs Frequency [THz]
    figure(5)
    
    a = plot(f./1E12,epsilon_r);
    hold on 
    b = plot(f./1E12,epsilon_i);
    plot(f./1E12,epsilon_r_mean);
    plot(f./1E12,epsilon_i_mean);
    hold off
    
    xlim([0.1 1]) 
    %ylim([-0.1 15]) 
    
    title_string = strcat('Dielectric constant SUBSTRATE (QHI,TETA)=', '(' ,string(QHI), ',', string(TETA),')');
    title(title_string);
    
    xlabel('Frequency[THz]')
    ylabel('\epsilon_r & \epsilon_i')
    legend([a,b],{'\epsilon_r','\epsilon_i'},'Location','southeast')
    
    saveas(gca, title_string);

%% Extimation with Brewster   
    %According to the measurements we made, we can consider a good
    %approximation for the Brewster angel the ange 73 degrees
    
    TETA_B = 73; %Degrees
    
    n_extimated = tan(TETA_B.*pi./180);
    
    n_approx = ones(L/2+1,1)*n_extimated;
    
    k_approx=1./(2.*n_approx).*(epsilon_i);
    k_mean_approx = k_approx(1:Cut_Index);
    k_mean_approx = ones(L/2+1,1).*mean(k_mean_approx);
    
%% Figure 6 - Index of Refraction vs Frequency [THz]
    figure(6)
    
    a = plot(f./1E12,n_approx);
    hold on 
    b = plot(f./1E12,k_approx);
    plot(f./1E12,k_mean_approx);
    hold off
    
    xlim([0.1 1]) 
    %ylim([-5 30]) 
    
    title_string = strcat('Index of refraction approximated (QHI,TETA)=', '(' ,string(QHI), ',', string(TETA),')');
    title(title_string);
    
    xlabel('Frequency[THz]')
    ylabel('n & k')
    legend([a,b],{'n','k'},'Location','southeast')    
    
    saveas(gca, title_string);
    
    
    
    
    
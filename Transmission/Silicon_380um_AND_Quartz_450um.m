%% Si 380um - Quartz 450um
% !n.b. The code save 8 graphs automatically when executed in the same
% folder
%% Presetting
clc;
clear all;  %clear the memmory
close all;  %close the previous plots and or data matrixes


%% Sample Selection - UnComment to choise ONE Sample

% Sample = 'Silicon';
Sample = 'Quartz';

%% Methodology Selection - UnComment to choise ONE Method

% Method = 'NoMethod';
Method = 'ZeroPadding';
% Method = 'BackgroundRemoval';

%% Sample Datas

switch Sample
    case 'Silicon'
        d1=380E-6; %m
        REF_path = 'REF_AIR_Silicon.txt';
        DATA_path = 'Silicon_380um.txt';
        TERALYZER_path = 'silicon_380um_D=380.csv';
    case 'Quartz'
        d1 = 450E-6; %m 
        REF_path = 'REF_AIR_Quartz.txt';
        DATA_path = 'Quartz_450um.txt';     
        TERALYZER_path = 'quartz_450um_D=450.csv';
    otherwise
        disp('ERROR: Select a valid Sample')
end
        

%% Datas 
%Thickness of the sample
                % d1=380E-6; %m   
                % %d1 = 450E-6; %m 
c=3E8; %m/s

%% Importing datas
%Reference
                % DATA1=importdata('REF_AIR_Silicon.txt'); %Reference
                % %DATA1=importdata('REF_AIR_Quartz.txt'); %Reference
                
DATA1=importdata(REF_path); %Reference

Time=DATA1(:,1);    %ps %Setting first column as time
Time=Time-Time(1);  %ps %Setting firt value to 0 and subsequent

SAMPLE=DATA1(:,2);  %a.u. %Setting second column as the amplitude
REF=DATA1(:,2);     %a.u. %Array for support in calculations

%Sample
                % DATA2=importdata('Silicon_380um.txt'); %Sample
                % %DATA2=importdata('Quartz_450um.txt'); %Sample
                
DATA2=importdata(DATA_path); %Sample

SUBSTRATE=DATA2(:,2);%a.u. %Setting second column as the amplitude



%% Uncomment this to apply Background Removal

switch Method
    case 'NoMethod'
    
    case 'ZeroPadding'
        TimeStep = Time(end) - Time(end-1);
        Time_index_S = 27; %ps
        Cut_Index_S = floor(Time_index_S./TimeStep);

        Sample_len = length(DATA1(:,2));
        SUBSTRATE = [DATA2(1:Cut_Index_S,2); zeros(Sample_len-(Cut_Index_S),1)];

            % I need the distance between the maximum of the SUBSTRATE and where i cut
            % Peak
        [M_S,num_S]=max(SUBSTRATE); %M_S = value of the max
                                    %num_s = index of M_S in SUBSTRATE                            
        Time_S=Time(num_S); %Time for the peak (M_S)
        [M_R,num_R]=max(REF); %M_S = value of the max
                                    %num_s = index of M_S in REF                            
        Time_R=Time(num_R); %Time for the peak (M_R)

        Time_index_R = Time_R + (Time_index_S-Time_S);
        Cut_Index_R = floor(Time_index_R./TimeStep); 

        REF = [DATA1(1:Cut_Index_R,2);zeros(Sample_len-(Cut_Index_R),1)];
    
    case 'BackgroundRemoval'
        TimeStep = Time(end) - Time(end-1);
        Time_index = 5; %ps %We are going to do a mean from 0 up to Time_index
        End_Index = floor(Time_index./TimeStep);

        Background_1 = mean(REF(1:End_Index));
        Background_2 = mean(SUBSTRATE(1:End_Index));

        REF = REF-Background_1;
        SUBSTRATE = SUBSTRATE-Background_2;
    otherwise disp('ERROR: insert a valid Method')
end


%% Figure 1 - E [a.u.] vs Time [ps]

    figure(1)
    plot(Time,REF) %Plot E of the reference vs Time (100ps)
    hold on
    plot(Time,SUBSTRATE) %Plot E of the reference vs Time (100ps)
    hold off
    
    title('THz Time Domain Spectra')
    xlabel('Time [ps]')
    ylabel('E Field [a.u]')
    legend({'AIR','SAMPLE'},'Location','southeast')
    
    %Save Figure    
    saveas(gca,'THz Time Domain Signal.png')

%% Variable change
    x=Time;
    y=REF;  
    z=SUBSTRATE;

%% FFT 
                          
    f_sampling = 1e12/(x(end)-x(end-1));    %1/Delta_t = Frequency Resolution
    L=length(y);    %Number of elements in the vector
    
    %f =(f_sampling*(0:(L/2))/L);
    %f=f';
    f = linspace(0,f_sampling/2,L/2+1)';    %Column of frequency axis

%% FFT Amplitudes    

    %FFT REFERENCE
    FFT_REF = fft(y,L)/L; % Normalized to L
    %AMPLITUDE OF FFT reference signal
    AFFT_REF = [abs(FFT_REF(1));2*abs(FFT_REF(2:L/2+1))];  %(Spalma la potenza) 
    %MAX AMPLITUDE OF FFT reference signal
    AFFT_REF_MAX=max(AFFT_REF);
    
    %FFT SAMPLE                      
    FFT_SUBSTRATE = fft(z,L)/L;
    %AMPLITUDE OF FFT sample signal
    AFFT_SUBSTRATE = [abs(FFT_SUBSTRATE(1));2*abs(FFT_SUBSTRATE(2:L/2+1))]; 
    %MAX AMPLITUDE OF FFT sample signal
    AFFT_SUBSTRATE_MAX=max(AFFT_SUBSTRATE);  
    
  
%% Figure 2 - Amplitude Spectra dB vs THz
    figure(2)
    
    semilogy(f./1E12,AFFT_REF)%Plotted in THz
    hold on
    semilogy(f./1E12,AFFT_SUBSTRATE)%Plotted in THz
    hold off
    
    xlim([0.1 6]) 
    title('Amplitude Spectra')
    xlabel('Frequency[THz]')
    ylabel('Amplitude [dB]')
    legend({'AIR','SAMPLE'},'Location','northeast')

    %Save Figure    
    saveas(gca,'Amplitude Spectra.png')
%% FFT Phases
    
    X2=FFT_REF;
    phase=atan2(imag(X2),real(X2))*180/pi;  %Four-quadrant arctan in degrees
    phase_single=phase(1:L/2+1);    %turns the plot into a one-sided plot
    
    %phase_single(2:end-1)=phase_single(2:end-1);
    
    %Unwrapping the phase
    Phase_REF= (180/pi)*unwrap(phase_single*pi/180);    %PHASE OF FFT in degrees

    
  
    X3=FFT_SUBSTRATE;
    phase=atan2(imag(X3),real(X3))*180/pi;
    phase_single=phase(1:L/2+1); %turns the plot into a one-sided plot
    
    %phase_single(2:end-1)=phase_single(2:end-1);

    %Unwrapping the phase
    Phase_SUBSTRATE= (180/pi)*unwrap(phase_single*pi/180); %PHASE OF FFT in degrees


%% Figure 3 - Phase [degrees] vs Frequency [THz]
    figure(3)
    
    plot(f./1E12,Phase_REF)
    hold on 
    plot(f./1E12,Phase_SUBSTRATE) 
    hold off
    
    xlim([0.1 2]) 
    title('FFT Phases')
    xlabel('Frequency[THz]')
    ylabel('Phase[Degree]')
    legend({'AIR','SAMPLE'},'Location','northeast')

    %Save Figure
    %saveas(gca,'FFT Phases.png')
    
    
%% EXTRACTING INDEX OF REFRACTION from given thickness
% % % d1=448E-6;%in m

w=2.*pi*f; %Vector of the angular frequencies

T=AFFT_SUBSTRATE./AFFT_REF;%TRANSFER FUNCTION FOR SUBSTRATE WITH RESPECT TO AIR
Delta_phase = Phase_REF-Phase_SUBSTRATE; %!n.b. This is in degree

d_min = d1 - 10E-6;
d_max = d1 + 10E-6;

%Real Part 
n=1+Delta_phase.*c./(180/pi)./(w.*d1); 

n_min = 1+Delta_phase.*c./(180/pi)./(w.*d_max); n_min(1) = n_min(2);
n_max = 1+Delta_phase.*c./(180/pi)./(w.*d_min); n_max(1) = n_max(2);
delta_n = (n_max-n_min)./2;

%Immaginary Part                                 
k=-1*(c./(w.*d1)).*log(T);%(extinction coefficient) 

%Better approx
%k=c./(w.*d1).*(log(4.*n./(n+1).^2)-log(abs(T)));

k_max = -1*(c./(w.*d_min)).*log(T); k_max(1) = k_max(2);
k_min = -1*(c./(w.*d_max)).*log(T); k_min(1) = k_min(2);
delta_k = (k_max-k_min)./2;


%% Figure 4 - Index of refractions vs Frequency [THz]
    figure(4)
    
    blue = [0.35 0.7 0.9];
    orange = [0.9,0.6,0];
    
    hold on
    
    %Filling the n graph
    fill([f'./1E12,fliplr(f'./1E12)], [n_min',fliplr(n_max')], blue, 'edgecolor', 'none', 'FaceAlpha', 0.5);
    %Filling the k graph
    fill([f./1E12;flipud(f./1E12)], [k_max;flipud(k_min)], orange, 'edgecolor', 'none', 'FaceAlpha', 0.5);
    
    p1 = plot(f./1E12,n, 'color', blue);
    p2 = plot(f./1E12,k, 'color', orange);
    hold off
    
    xlim([0.1 2]) 
    ylim([-0.1 5]) 
    title('Index of refraction SUBSTRATE')
    xlabel('Frequency[THz]')
    ylabel('n&k')
    legend([p1,p2],{'n','k'},'Location','northeast')
    
    %Save Figure
    saveas(gca,'Index of refraction SUBSTRATE.png')

    
%% Extraction of other parameters:


%% Complex Dielectric Constant

    %Real Part 
    epsilon_r = n.^2-k.^2;
    
    epsilon_r_max = n_max.^2-k_max.^2;
    epsilon_r_min = n_min.^2-k_min.^2;
    
    %Immaginary Part                                 
    epsilon_i = 2*n.*k;
    
    epsilon_i_max = 2*n_max.*k_max;
    epsilon_i_min = 2*n_min.*k_min;
    
%% Figure 5 - Dielectric constant vs Frequency [THz]
    figure(5)
    
    hold on 
    
    fill([f'./1E12,fliplr(f'./1E12)], [epsilon_r_max', fliplr(epsilon_r_min')], blue, 'edgecolor', 'none', 'FaceAlpha', 0.5);
    fill([f'./1E12,fliplr(f'./1E12)], [epsilon_i_max', fliplr(epsilon_i_min')], orange, 'edgecolor', 'none', 'FaceAlpha', 0.5);
    
    p3 = plot(f./1E12,epsilon_r, 'color', blue);
    p4 = plot(f./1E12,epsilon_i, 'color', orange);
    hold off
    xlim([0.1 2]) 
    ylim([-0.1 15]) 
    title('Dielectric constant SUBSTRATE')
    xlabel('Frequency[THz]')
    ylabel('\epsilon_r & \epsilon_i')
    legend([p3,p4], {'\epsilon_r','\epsilon_i'},'Location','northeast')

    %Save Figure    
    saveas(gca,'Dielectric Constant.png')    
    
%% Complex Conducibility

    epsilon_0 = 8.85E-12; %C^2/(Nm^2)
    epsilon_infty = 4;
    %Real Part 
    sigma_r = (epsilon_0.*w.*epsilon_i)./100; %Siemens/cm
    %Immaginary Part                                 
    %sigma_i = epsilon_0.*w.*(epsilon_infty-epsilon_r)./100;
    
 %% Figure 6 - Conducibility vs Frequency [THz]
    figure(6)
    
    plot(f./1E12,sigma_r)
%     hold on 
%     plot(f./1E12,sigma_i)
%     hold off
    xlim([0.1 2]) 
    ylim([-0.1 5]) 
    title('Complex conducibility SUBSTRATE')
    xlabel('Frequency[THz]')
    ylabel('\sigma_r [\Omega^{-1}cm^{-1}]')
    legend({'\sigma_r'},'Location','northeast') 

    %Save Figure
    saveas(gca,'Complex conducibility SUBSTRATE.png')
    
    
% %% Absorbance coefficient
%     
% %     R = abs(((n+j.*k-1)./(n+j.*k+1)).^2);
% %     Transmittance = 1 - R;
% 
% %     alpha = ((-2/d1).*log(T))./100;
%   alpha = (2.*w.*k)./c./100;
%     
%  %% Figure 7 - Absorbance coefficient vs Frequency [THz]
%     figure(7)
%     
%     plot(f./1E12,alpha)
%     
%     xlim([0.1 2]) 
%     %ylim([-0.1 5]) 
%     title('Absorbance coefficient SUBSTRATE')
%     xlabel('Frequency[THz]')
%     ylabel('\alpha')
%     legend({'\alpha'},'Location','southeast')
%     
%     %Save Figure    
%     exportgraphics(gca,'Absorbance coefficient SUBSTRATE.png','Resolution',300)
    
    
%% Comparison with Teralyzer    

%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: C:\Users\leona\Desktop\THz_transmission_data_onSilicon and quartz_with MATLAB_Code\silicon_380um_D=380.csv
%
% Auto-generated by MATLAB on 03-Apr-2021 12:36:23

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 29);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["freq", "ref_ind", "kappa", "alpha", "epsilon_r", "epsilon_i", "ref_ind_S", "kappa_S", "alpha_S", "alpha_max", "epsilonS_r", "epsilonS_i", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29"];
opts.SelectedVariableNames = ["freq", "ref_ind", "kappa", "alpha", "epsilon_r", "epsilon_i", "ref_ind_S", "kappa_S", "alpha_S", "alpha_max", "epsilonS_r", "epsilonS_i"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29"], "EmptyFieldRule", "auto");

% Import the data
TERALYZER = readtable(TERALYZER_path, opts);

%% Convert to output type
TERALYZER = table2array(TERALYZER);

%% Clear temporary variables
clear opts

Frequencies = TERALYZER(:,1)/(1E12);
Real_n = TERALYZER(:,2);
Imaginary_n = TERALYZER(:,3);
Absorbance = TERALYZER(:,4);
Real_epsilon = TERALYZER(:,5);
Imaginary_epsilon = TERALYZER(:,6);

%% Plot to compare datas - Index of Refraction
figure(8)

    hold on
    plot(Frequencies,Real_n)
    plot(Frequencies,Imaginary_n)
    
    plot(f./1E12,n)
    plot(f./1E12,k)
    hold off
    
    xlim([0.3 1]) 
    ylim([-0.1 5]) 

    title('Index of refraction - Comparison')
    xlabel('Frequency[THz]')
    ylabel('n & k')
    legend({'n_{Teralyzer}', 'k_{Teralyzer}', 'n', 'k'}, 'Location', 'northeast')
    
    %Save Figure    
    saveas(gca,'Index of refraction - Comparison.png')
    
figure(9)

    hold on
    plot(Frequencies,Real_epsilon)
    plot(Frequencies,Imaginary_epsilon)
    
    plot(f./1E12,epsilon_r)
    plot(f./1E12,epsilon_i)
    hold off
   
    xlim([0.3 1]) 
    ylim([-0.1 13]) 

    title('Dielectric Constant - Comparison')
    xlabel('Frequency[THz]')
    ylabel('\epsilon_r & \epsilon_i')
    legend({'\epsilon_r_{Teralyzer}', '\epsilon_i_{Teralyzer}', '\epsilon_r', '\epsilon_i'}, 'Location', 'east') 
    
    %Save Figure
    saveas(gca,'Dielectric Constant - Comparison.png')
    
% figure(10)
%     hold on
%     plot(Frequencies,Absorbance)
%     plot(f./1E12,alpha)
%     hold off
%    
%     xlim([0.3 1]) 
%     ylim([-0.1 13]) 
%     
%     %Save Figure
%     exportgraphics(gca,'fig10.png','Resolution',300)
%% Saving Plots

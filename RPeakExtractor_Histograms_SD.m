%% R Peak Detection, Synthesized (for all subjects of interest)
close all
%EDIT
%load subjects; note - set path to file location of subjects
%load HEPv1_020_20200213_040427.mat %subject 20
%load HEPv1_021_20200220_024026.mat %subject 21
load HEPv1_022_20200224_031331.mat %subject 22
%EDIT

%matrix of ECG variable names per subject
all_data_ecg = {HEPv1_022_20200224_031331mff1ECG
HEPv1_022_20200224_031331mff2ECG
HEPv1_022_20200224_031331mff3ECG
HEPv1_022_20200224_031331mff4ECG
HEPv1_022_20200224_031331mff5ECG
HEPv1_022_20200224_031331mff6ECG
HEPv1_022_20200224_031331mff7ECG
HEPv1_022_20200224_031331mff8ECG
HEPv1_022_20200224_031331mff9ECG
HEPv1_022_20200224_031331mff10ECG
HEPv1_022_20200224_031331mff11ECG
HEPv1_022_20200224_031331mff12ECG
HEPv1_022_20200224_031331mff13ECG
HEPv1_022_20200224_031331mff14ECG
HEPv1_022_20200224_031331mff15ECG
HEPv1_022_20200224_031331mff16ECG
HEPv1_022_20200224_031331mff17ECG
HEPv1_022_20200224_031331mff18ECG
HEPv1_022_20200224_031331mff19ECG
HEPv1_022_20200224_031331mff20ECG
HEPv1_022_20200224_031331mff21ECG
HEPv1_022_20200224_031331mff22ECG
HEPv1_022_20200224_031331mff23ECG
HEPv1_022_20200224_031331mff24ECG
HEPv1_022_20200224_031331mff25ECG
HEPv1_022_20200224_031331mff26ECG
HEPv1_022_20200224_031331mff27ECG
HEPv1_022_20200224_031331mff28ECG
HEPv1_022_20200224_031331mff29ECG};

%% Histogram Generation
figure
%% CONDITION 1
subplot(2,2,1)
hold on
R_R = [];
HR = 0;
Time_tot = 0;
for i = 1:7
    ECGmat = cell2mat(all_data_ecg(i)); time_s = make_time_s(ECGmat,EEGSamplingRate); %generate time scale
    R_peaks = find_Rpeaks(time_s,ECGmat); %<<---
    R_R = [R_R, diff(time_s(R_peaks))];
    HR = length(R_peaks) + HR; %calculate beats
    Time_tot = Time_tot + time_s(end); %total time elasped
end
HR = HR/Time_tot*60; %calculate heart rate in bpm
fprintf("\nC1 \nHR (bpm): %f",HR)
histogram(R_R.*1000); %in ms
fitdist(R_R'.*1000,'Normal')
xlabel('R-R in ms');
title('Condition 1')
RMSSD = rms(diff(R_R).*1000); %in ms 
fprintf("RMSSD (ms): %f",RMSSD);

hold off
%% CONDITION 2
subplot(2,2,2)
hold on
R_R = [];
HR = 0;
Time_tot = 0;
for i = 8:14
    ECGmat = cell2mat(all_data_ecg(i)); time_s = make_time_s(ECGmat,EEGSamplingRate); %generate time scale
    R_peaks = find_Rpeaks(time_s,ECGmat); %<<---
    R_R = [R_R, diff(time_s(R_peaks))];
    HR = length(R_peaks) + HR; %calculate beats
    Time_tot = Time_tot + time_s(end); %total time elasped
end
HR = HR/Time_tot*60; %calculate heart rate in bpm
fprintf("\nC2 \nHR (bpm): %f",HR)
histogram(R_R.*1000); %in ms
fitdist(R_R'.*1000,'Normal')
xlabel('R-R in ms');
title('Condition 2')
RMSSD = rms(diff(R_R).*1000); %in ms 
fprintf("RMSSD (ms): %f",RMSSD);

hold off
%% CONDITION 3
subplot(2,2,3)
hold on
R_R = [];
HR = 0;
Time_tot = 0;
for i = 16:22
    ECGmat = cell2mat(all_data_ecg(i)); time_s = make_time_s(ECGmat,EEGSamplingRate); %generate time scale
    R_peaks = find_Rpeaks(time_s,ECGmat); %<<---
    R_R = [R_R, diff(time_s(R_peaks))];
    HR = length(R_peaks) + HR; %calculate beats
    Time_tot = Time_tot + time_s(end); %total time elasped
end
HR = HR/Time_tot*60; %calculate heart rate in bpm
fprintf("\nC3 \nHR (bpm): %f",HR)
histogram(R_R.*1000); %in ms
fitdist(R_R'.*1000,'Normal')
xlabel('R-R in ms');
title('Condition 3')
RMSSD = rms(diff(R_R).*1000); %in ms 
fprintf("RMSSD (ms): %f",RMSSD);

hold off
%% CONDITION 4
subplot(2,2,4)
hold on
R_R = [];
HR = 0;
Time_tot = 0;
for i = 23:length(all_data_ecg)
    ECGmat = cell2mat(all_data_ecg(i)); time_s = make_time_s(ECGmat,EEGSamplingRate); %generate time scale
    R_peaks = find_Rpeaks(time_s,ECGmat); %<<---
    R_R = [R_R, diff(time_s(R_peaks))];
    HR = length(R_peaks) + HR; %calculate beats
    Time_tot = Time_tot + time_s(end); %total time elasped
end
HR = HR/Time_tot*60; %calculate heart rate in bpm
fprintf("\nC4 \nHR (bpm): %f",HR)
histogram(R_R.*1000); %in ms
fitdist(R_R'.*1000,'Normal')
xlabel('R-R in ms');
title('Condition 4')
RMSSD = rms(diff(R_R).*1000); %in ms 
fprintf("RMSSD (ms): %f",RMSSD);

hold off


%% Helper Functions:

%% Generate the time scale
function time_s = make_time_s(ECGmat, EEGSamplingRate)
   
    EEG_fs = EEGSamplingRate; %sampling frequency equal to EEG rate
    
    % Where sample_n = length(HEP)
    sample_n = length(ECGmat);
    
    EEG_ts = 1/EEG_fs; %in s
    time_s = linspace(0,EEG_ts*(sample_n-1),sample_n); %time scale in s
    
end

% %%
% function y = process_signal_noisy(x)
% end

%% Find R Peaks
function R_index = find_Rpeaks(tm,ecgsig)
        %% Source info: https://www.mathworks.com/help/wavelet/ug/r-wave-detection-in-the-ecg.html
    % References
        % 
        % Goldberger A. L., L. A. N. Amaral, L. Glass, J. M. Hausdorff, P. Ch. Ivanov, R. G. Mark, 
        % J. E. Mietus, G. B. Moody, C-K Peng, H. E. Stanley. "PhysioBank, PhysioToolkit, and PhysioNet: 
        % Components of a New Research Resource for Complex Physiologic Signals." Circulation 101. Vol.23, 
        % e215-e220, 2000. http://circ.ahajournals.org/cgi/content/full/101/23/e215
        % 
        % Moody, G. B. "Evaluating ECG Analyzers". http://www.physionet.org/physiotools/wfdb/doc/wag-src/eval0.tex
        % 
        % Moody G. B., R. G. Mark. "The impact of the MIT-BIH Arrhythmia Database." IEEE Eng in Med and Biol. 
        % Vol. 20, Number 3, 2001), pp. 45-50 .
    
    %% Load and plot an ECG waveform where the R peaks of the QRS complex have
    % been annotated by two or more cardiologists.
        
    %     load mit200 %(data)
    %     figure
    %     plot(tm,ecgsig) %EX; tm = mme_s; ecgsig = dataECG;
    %     hold on
    %     plot(tm(ann),ecgsig(ann),'ro')
    %     xlabel('Seconds')
    %     ylabel('Amplitude')
    %     title('Subject - MIT-BIH 200') %EX; title('Subject - 20 - Community, Epoch 1')
    
    
%     %load mit200
%     figure
%     plot(tm,ecgsig)
%     hold on
%     plot(tm,ecgsig)
%     xlabel('Seconds')
%     ylabel('Amplitude')
%     title('Subject - 20 - Community, Matrix 2')
    
    %% DNC - The 'sym4' wavelet resembles the QRS complex, which makes it a good choice 
    % for QRS detection. 
    % To illustrate this more clearly, extract a QRS complex and plot the result 
    % with a dilated and translated 'sym4' wavelet for comparison.
    
        % qrsEx = ecgsig(4560:4810);
        % [mpdict,~,~,longs] = wmpdictionary(numel(qrsEx),'lstcpt',{{'sym4',3}});
        % figure
        % plot(qrsEx)
        % hold on
        % plot(2*circshift(mpdict(:,11),[-2 0]),'r')
        % axis tight
        % legend('QRS Complex','Sym4 Wavelet')
        % title('Comparison of Sym4 Wavelet and QRS Complex')
    
    %% DNC - Use the maximal overlap discrete wavelet transform (MODWT) to enhance the R peaks 
    % in the ECG waveform. The MODWT is an undecimated wavelet transform, 
    % which handles arbitrary sample sizes.
    
    % First, decompose the ECG waveform down to level 5 using the default 'sym4' wavelet.
    % Then, reconstruct a frequency-localized version of the ECG waveform using only the 
    % wavelet coefficients at scales 4 and 5. The scales correspond to the following approximate
    % frequency bands.
    % 
    % Scale 4 -- [11.25, 22.5) Hz
    % 
    % Scale 5 -- [5.625, 11.25) Hz.
    % 
    % This covers the passband shown to maximize QRS energy.
    
        % wt = modwt(ecgsig,5);
        % wtrec = zeros(size(wt));
        % wtrec(4:5,:) = wt(4:5,:);
        % y = imodwt(wtrec,'sym4');
    
    wt = modwt(ecgsig,5);
    wtrec = zeros(size(wt));
    wtrec(4:5,:) = wt(4:5,:);
    y = imodwt(wtrec,'sym4');
    
    %% Use the squared absolute values of the signal approximation built from the
    % wavelet coefficients and employ a peak finding algorithm to identify the R peaks.
    
    % If you have Signal Processing Toolbox™, you can use findpeaks to locate the peaks.
    % Plot the R-peak waveform obtained with the wavelet transform annotated with the 
    % automatically-detected peak locations.
        
        % y = abs(y).^2;
        % [qrspeaks,locs] = findpeaks(y,tm,'MinPeakHeight',0.35,...
        %     'MinPeakDistance',0.150);
        % figure
        % plot(tm,y)
        % hold on
        % plot(locs,qrspeaks,'ro')
        % xlabel('Seconds')
        % title('R Peaks Localized by Wavelet Transform with Automatic Annotations')
    
    y = abs(y).^2;
%     %% Debug
%     figure
%     hold on
%     plot(y)
%     hold off
%     %pause%
%     %%
    
    [qrspeaks,locs] = findpeaks(y,tm,'MinPeakHeight',1.85E5,...
        'MinPeakDistance',0.150); %minpeakheight 2.25E5
%     figure
%     plot(tm,y)
%     hold on
%     plot(locs,qrspeaks,'ro')
%     xlabel('Seconds')
%     title('R Peaks Localized by Wavelet Transform with Automatic Annotations')
%     hold off
    %pause%
    %% Add the expert annotations to the R-peak waveform. Automatic peak detection times are considered 
    % accurate if within 150 msec of the true peak (±75 msec)
        % plot(tm(ann),y(ann),'k*')
        % title('R peaks Localized by Wavelet Transform with Expert Annotations')
%     %% DNC - Automated Peak Detection Compared with Raw Signal 
%     figure
%     plot(tm,ecgsig)
%     hold on
%     %edit%
%     locs_y = sqrt(qrspeaks); %for depiction purposes
%     %edit%
%     plot(locs,locs_y,'ro') % edit: plot(locs,qrspeaks,'ro')
%     xlabel('Seconds')
%     ylabel('Amplitude')
%     title('R Peaks Localized by Automatic Annotations')
%     hold off
    
    %% Convert Find peaks into indices
    R_index = ones(1,length(locs)).*-1;
    for i = 1:length(locs)
        R_index(i) = find(tm==locs(i)); %Find index
    end
%     %% Create a Table s.t. 
%     % R_index is "R Peak Indices," locs is "R Peak Times (s)," where time is calculated by 
%     % tm = time_s = linspace(0,EEG_ts*sample_n,sample_n); %time scale in s
%     varNames = {'R Peak Indices','R Peak Times (s)'}; %,'Time Vector'};
%     R_index = R_index';
%     locs = locs';
%     T = table(R_index,locs,'VariableNames',varNames);
%     
%     varNames = {'Time Vector'};
%     tm = tm';
%     T2 = table(tm,'VariableNames',varNames);
%     %%edit
%     filename = 'R_index_Mat2_Sub20.txt';
%     %%edit
%     writetable(T,filename);
%     
%     
%     %%edit
%     filename = 'TimeVectorUsed_Mat2_Sub20.txt';
%     %%edit
%     writetable(T2,filename);
end

%%
% function y = find_Twave(tm) 
% end
% %%
% function y = process_for_t(x)
% end

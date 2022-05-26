%% R Peak Detection, Synthesized (for all subjects of interest)
close all

f1 = figure("Name",'Distribution of Peaks');
f2 = figure("Name",'Peak Location vs Cleaned');
%EDIT
%load subjects; note - set path to file location of subjects
%load HEPv1_020_20200213_040427.mat %subject 20
load HEPv1_021_20200220_024026.mat %subject 21
%load HEPv1_022_20200224_031331.mat %subject 22
%EDIT

%Targeted Matrices
%Instruction: Copy variable names from workspace into excel. 
%Filter for ECG values in COLUMN. 
%May need to generate txt/csv. Paste between brackets below. 
%Use all ECGs!
all_data_ecg = {HEPv1_021_20200220_024026mff1ECG
HEPv1_021_20200220_024026mff2ECG
HEPv1_021_20200220_024026mff3ECG
HEPv1_021_20200220_024026mff4ECG
HEPv1_021_20200220_024026mff5ECG
HEPv1_021_20200220_024026mff6ECG
HEPv1_021_20200220_024026mff7ECG
HEPv1_021_20200220_024026mff8ECG
HEPv1_021_20200220_024026mff9ECG
HEPv1_021_20200220_024026mff10ECG
HEPv1_021_20200220_024026mff11ECG
HEPv1_021_20200220_024026mff12ECG
HEPv1_021_20200220_024026mff13ECG
HEPv1_021_20200220_024026mff14ECG
HEPv1_021_20200220_024026mff15ECG
HEPv1_021_20200220_024026mff16ECG
HEPv1_021_20200220_024026mff17ECG
HEPv1_021_20200220_024026mff18ECG
HEPv1_021_20200220_024026mff19ECG
HEPv1_021_20200220_024026mff20ECG
HEPv1_021_20200220_024026mff21ECG
HEPv1_021_20200220_024026mff22ECG
HEPv1_021_20200220_024026mff23ECG
HEPv1_021_20200220_024026mff24ECG
HEPv1_021_20200220_024026mff25ECG
HEPv1_021_20200220_024026mff26ECG
HEPv1_021_20200220_024026mff27ECG
HEPv1_021_20200220_024026mff28ECG
HEPv1_021_20200220_024026mff29ECG};

%% Histogram Generation
%% CONDITION 1: Find R_R of all peaks
%if
first = 1; %matrices needed
last = 7;
R_R = [];
HR = 0;
Time_tot = 0;
R_peaks = [];

%Trial Indices:
End_Trial = []; %ending indices of each trial in R_peaks_clean
for i = first:last
    ECGmat = cell2mat(all_data_ecg(i)); time_s = make_time_s(ECGmat,EEGSamplingRate); %generate time scale
    R_peaks_trial = find_Rpeaks(time_s,ECGmat); %<<---indices per trial
    R_peaks = [R_peaks R_peaks_trial]; %<<---indices total condition
    R_R = [R_R, diff(time_s(R_peaks_trial))];
    HR = length(R_peaks_trial) + HR; %calculate beats
    Time_tot = Time_tot + time_s(end); %total time elasped
   
    End_Trial = [End_Trial length(R_peaks_trial)];

    %%check or plot the third matrix *******************************
    if i == 3 
        hold on
        figure(f2)     
        y = zeros(1,length(R_peaks_trial));
        plot(R_peaks_trial,y,'ro','Marker', 'o','Color','k'); %Raw peaks
        %Raw
        plot(ECGmat); %Raw peaks
        hold off
    end
    %%check or plot the third matrix *******************************

end

%% Statistical Summary
hold on
figure(f1)

HR = HR/Time_tot*60; %calculate heart rate in bpm
fprintf("\nC1 \nHR (bpm): %f",HR)
histogram(R_R.*1000); %in ms
fitdist(R_R'.*1000,'Normal')
xlabel('R-R in ms');
title('Condition 1')
RMSSD = rms(diff(R_R).*1000); %in ms 
fprintf("RMSSD (ms): %f",RMSSD);

hold off

%% Clean R_R peaks
%Indices of outliers
outlier = 3*std(R_R); %R_R index, backtracked to the R matrix; uses 99.7 %CI
mu = mean(R_R);

FlR = find(R_R>(mu+outlier)); %point to the right
FlR = FlR+1; %convert to R locations; peak locations in Rpeaks that are outliers 
FlL = find(R_R<(mu-outlier)); %point to the left (as is); peak locations in Rpeaks that are outliers

Fl = [FlL FlR];
Fl = sort(Fl); %indices of outliers 

count=1;
count2=1;
R_peaks_clean = []; %R_peaks of a full condition without outliers 
for i = 1:length(R_peaks)
    if ~isempty(Fl)&&count <= length(Fl)&&(i==Fl(count))
        count = count + 1; %skip assignment
    else
        R_peaks_clean(count2) = R_peaks(i); 
        count2 = count2 + 1;
        %skip assignment
    end
end

%% Print all peaks
% T = table(R_peaks_clean);
filename = sprintf("R_Peak_Clean_C%d_AllTrials.txt",1); 
% writetable(T,filename);

fileID = fopen(filename,'w');
%fprintf(fileID,'%6s %12s\n','x','exp(x)'); %title
fprintf(fileID,'%d\n',R_peaks_clean);
fclose(fileID);

%% Print by matrix/trial/mff
loc=1;%track of matrix
dat=[];
for i = 1:length(End_Trial)
    if ~(i==1)
        dat = R_peaks_clean(loc:End_Trial(i)+loc-1);
        loc = loc + End_Trial(i)
    else
        dat = R_peaks_clean(loc:End_Trial(i));
        loc = End_Trial(i) + 1
    end

    %%check or plot the third matrix *******************************
    if i==3   
        hold on
        figure(f2)     
        y = zeros(1,length(dat));
        plot(dat,y,'ro','Marker','*','Color','b');
        hold off
    end
    %%optional plot *******************************


    %Create Text File for Matrix i

    %T = table(R_peaks_clean);
    filename = sprintf('R_Peak_Clean_C%d_Trial%d.txt',1,i); 
    %writetable(T,filename);

    fileID = fopen(filename,'w');
    %fprintf(fileID,'%6s %12s\n','x','exp(x)'); %title
    fprintf(fileID,'%d\n',dat);
    fclose(fileID);
end

%% Labelling Trials
%% CONDITION 2
%if else
%% CONDITION 3
%if else
%% CONDITION 4
%else
%subplot(2,2,4)

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

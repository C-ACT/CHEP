

%% Community Sample Script 
% Date of Creation: 02-08-2022 - FINAL; Edited 06-14-2022
% Author: Samantha Pitts @ Mount Sinai West Waters Lab
%Additions: Visualizing performance of algorithm
% References: 
        % Source info: https://www.mathworks.com/help/wavelet/ug/r-wave-detection-in-the-ecg.html
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

 % Conditions:
 % Community sample trials:
    % 
    % Resting eyes open: 7 trials 
    % Resting eyes closed: 7 trials 
    % 
    % HBP Eyes Open: 7 trials; 
    % HBP Eyes Closed: 7 trials

    %% Command Window Signal Quality Check
    %T = timetable(Var','SampleRate',Fs) where Fs = 1000 in NetStation
    %% Init.

Total_Trials = length( ...
   all_data_ecg);
Trial_Lengths = []; %size of each trial for timekeeping
R_peaks_trial_start = [1];%keep track of trials from list of peaks
End_Trial = []; %ending indices of each trial in R_peaks_clean
Trial_Label = []; %
Epoch_No = [];
T = [];


%% R Peak Detection, Synthesized (for all subjects of interest)
close all


%% EDIT
%load subjects; note - set path to file location of subjects
%cd '/Users/samanthapitts/Desktop/Research | Mount Sinai/[Project] - HEP/Matlab File/Community Sample HEP/Data/'
%load Files/HEPv1_020_20200213_040427.mat %subject 20
%load HEPv1_021_20200220_024026.mat %subject 21
%load HEPv1_022_20200224_031331.mat %subject 22
[file, path] = uigetfile;
cd(path)
load(file)
%EDIT

% %Select target folder
% [~, path] = uigetfile;
% cd(path)
%% Targeted Matrices
%Instruction: Copy variable names from workspace into excel. 
%Filter for ECG values in COLUMN. Use:  Variable names - CHEP 1.csv
%Convert names into "Plain Text." Paste between brackets below. 
%Use all ECGs! (1-29)
all_data_ecg = {HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359
HEPv1_015_20191210_022359};


%%EDIT
Subject_No = 15;
Condition_Matrix_Start = 2:7:Total_Trials; %C1 Starts @ Mat 2, C2 Starts 7 after,... until last Mat 29
plots = 28; %28 unique trials for peak detection
f1 = figure('Name','R Peaks Localized by Abs Square Data');
%f2 = figure('Name','R Peaks Localized by Abs Square Data - with Raw Pairing');

%Method: Fourth Order Butterworth bp filter to isolate reasonable HR. DWT Method.
f2 = figure('Name','R Peaks Localized by Processed Peak Prominence');
%%EDIT

%% Find R_R of all peaks

%Throw-away trial for NetStation Timing
Trial_Lengths = [length(cell2mat(all_data_ecg(1)))]; %1st matrix place holder

% PER CONDITION "A"
Num_Conditions = length(Condition_Matrix_Start); %Community has 4 conditions recorded
R_peaks = [];
p = 0; %plot number
for j = 1:Num_Conditions
    
    first = Condition_Matrix_Start(j); %matrix number per condition
    if (j==Num_Conditions) 
        last = Total_Trials; %last matrix belonging to condition
    else
        last = Condition_Matrix_Start(j+1)-1;
    end

    % PER TRIAL IN Condition "A"
    for i = first:last %for each trial
        p = p+1;%plot number
        ECGmat = cell2mat(all_data_ecg(i)); 
        time_s = make_time_s(ECGmat,EEGSamplingRate); %generate time scale
        R_peaks_trial = find_Rpeaks(time_s,ECGmat,plots,f1,f2,p); %<<---indices per trial
        
        R_peaks = [R_peaks R_peaks_trial]; %<<---indices all 

        if ~(j==4 && i==7) %(last mat) 
            R_peaks_trial_start = [R_peaks_trial_start, length(R_peaks)+1];
        end
        
        Trial_Lengths = [Trial_Lengths length(ECGmat)]; %s.t. [applytogroup1(length mat 1), applytogroup2+gr1...]
        
        x = "error";
        if j == 1 %per condition, for each trial
             x = "MWEO";
        elseif j == 2    
             x = "MWEC";
        elseif j == 3
             x = "HBEO";
        else
             x = "HBEC";
        end
        z = repmat(x,1,length(R_peaks_trial));
        Trial_Label = [Trial_Label, z];

        z2 = repmat(i,1,length(R_peaks_trial));
        Epoch_No = [Epoch_No, z2];

    end
end
%% Format Peak Via NetStation
R_peaks_NS = R_peaks;

%% ADDENDUM: Time for NS
R_peaks_NS; %locations of R peaks;
R_peaks_NS = (R_peaks_NS - 1); %NS adjusted peaks (for visual purposes);
time_NS = seconds(R_peaks_NS/1000); %into sec per Sample Rate = 1000Hz
time_NS_formatted = duration(time_NS,'Format','hh:mm:ss.SSS'); %formatted NS times

%% Print all peaks
% T = table(R_peaks_clean);
filename = sprintf("NS_R_Peaks_Subject%d_All.txt",Subject_No); 
% writetable(T,filename);

fileID = fopen(filename,'w');
%fprintf(fileID,'%6s %12s\n','x','exp(x)'); %title
fprintf(fileID,'%d\n',R_peaks_NS);
fclose(fileID);
%% Print NetStation Times
filename = sprintf("NSEpoch_String_R_Peaks_Subject%d_All.txt",Subject_No); 

fileID = fopen(filename,'w');
time = string(time_NS_formatted);
for i = 1:length(Epoch_No)
    fprintf(fileID,'_[%d] %s\n',Epoch_No(i),time(i));
end
fclose(fileID);


%% Print all peaks
filename = sprintf("NSEpoch_R_Peaks_Subject%d_All.txt",Subject_No); 

fileID = fopen(filename,'w');
fprintf(fileID,'%d\n',R_peaks_NS);
fclose(fileID);

%% Print Labels
filename = sprintf("NSEpoch_Label_Subject%d_All.txt",Subject_No); 

fileID = fopen(filename,'w');
fprintf(fileID,'%s\n',Trial_Label);
fclose(fileID);

%% save figure
filename = sprintf("Sub%d_Peak_by_epoch",Subject_No); 
savefig(f2,filename);

%% Helper Functions:
%% Generate the time scale
function time_s = make_time_s(ECGmat, EEGSamplingRate)
   
    EEG_fs = EEGSamplingRate; %sampling frequency equal to EEG rate
    
    % Where sample_n = length(HEP)
    sample_n = length(ECGmat);
    
    EEG_ts = 1/EEG_fs; %in s
    time_s = linspace(0,EEG_ts*(sample_n-1),sample_n); %time scale in s
    
end
%% Find R Peaks
function R_index = find_Rpeaks(tm,ecgsig,plots,f1,f2,p)
%finding the peaks of n trials is equivalent to n plots (Ex:28)
%establish two figures f1, f2 to save plot analysis
    %% Load and plot an ECG waveform where the R peaks of the QRS complex have
    % been annotated by two or more cardiologists.
    
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
    
    %% Preprocessing 
    load EKG_Butterworth_BP.mat
    ecgsig = filter(EKG_Butterworth_BP,ecgsig);

    %% DNC - Use the maximal overlap discrete wavelet transform (MODWT) to enhance the R peaks 
    % in the ECG waveform. The MODWT is an undecimated wavelet transform, 
    % which handles arbitrary sample sizes.
    wt = modwt(ecgsig,5);
    wtrec = zeros(size(wt));
    wtrec(4:5,:) = wt(4:5,:);
    y = imodwt(wtrec,'sym4');
    
    %% Use the squared absolute values of the signal approximation built from the
    % wavelet coefficients and employ a peak finding algorithm to identify the R peaks.
    
    % If you have Signal Processing Toolbox™, you can use findpeaks to locate the peaks.
    % Plot the R-peak waveform obtained with the wavelet transform annotated with the 
    % automatically-detected peak locations.
    

    %% Prior 
    %y = abs(y).^2;
    %% Parameters  
    [qrspeaks,locs] = findpeaks(y,tm,'MinPeakHeight',1,...
        'MinPeakDistance',0.150,'MinPeakProminence',500); %minpeakheight 2.25E5 663+500 =1000
    
    %% For Summary Plots in single file
    n = floor(plots/7);
    check = mod(n,7);
    if check~=0
        n = n+1;
    end

    %% plots (first one for abs squared peak detection)
%     figure(f1)
%     subplot(7,n,p)
%     hold on
%     plot(tm,y)
%     plot(locs,qrspeaks,'ro')
%     xlabel('Seconds')
%     %title('R Peaks Localized by Wavelet Transform with Automatic Annotations')
%     hold off

    figure(f2)
    subplot(7,4,p)
    hold on
    plot(tm,ecgsig)
    idx=[];
    for k= 1:length(locs)
        idx(k) = find(tm==locs(k));
    end
    plot(locs,ecgsig(idx),'ro')
    xlabel('Seconds')
    %title('R Peaks Localized by Raw Data')
    hold off    
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

end

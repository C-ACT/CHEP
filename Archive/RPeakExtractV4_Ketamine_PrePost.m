%% Ketamine Sample Script with TIME
% Date of Creation: 03-02-2022; Last update: 03/03/2022
% Author: Samantha Pitts @ Mount Sinai West Waters Lab
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
 % Ketamine trials File PRE:
    % Blank Matrix: 1 trial/matrix
    % MMN: 1 trials (2)
    % ReO: 1 trials (3)
    % ReCl: 1 trials (4)
    % HBP Eyes Open: 7 trials (5-11)

%% R Peak Detection, Synthesized (for all subjects of interest)
close all

% f1 = figure("Name",'Distribution of Peaks');
% f2 = figure("Name",'Peak Location vs Cleaned');

%EDIT
%load subjects; note - set path to file location of subjects
% load PRE_MMN_HEP_REST_01_06_22_20220106_013811.mat %subject
%EDIT

%Targeted Matrices
%Instruction: Copy variable names from workspace into excel. 
%Filter for ECG values in COLUMN. 
%May need to generate txt/csv. Paste between brackets below. 
%Use all ECGs!
all_data_ecg = {KPILOT_002_05_POST_20220126_032514mff1ECG
KPILOT_002_05_POST_20220126_032514mff2ECG
KPILOT_002_05_POST_20220126_032514mff3ECG
KPILOT_002_05_POST_20220126_032514mff4ECG
KPILOT_002_05_POST_20220126_032514mff5ECG
KPILOT_002_05_POST_20220126_032514mff6ECG
KPILOT_002_05_POST_20220126_032514mff7ECG
KPILOT_002_05_POST_20220126_032514mff8ECG
KPILOT_002_05_POST_20220126_032514mff9ECG
KPILOT_002_05_POST_20220126_032514mff10ECG
KPILOT_002_05_POST_20220126_032514mff11ECG};
% {PRE_MMN_HEP_REST_01_06_22_20220106_013811mff1ECG
% PRE_MMN_HEP_REST_01_06_22_20220106_013811mff2ECG
% PRE_MMN_HEP_REST_01_06_22_20220106_013811mff3ECG
% PRE_MMN_HEP_REST_01_06_22_20220106_013811mff4ECG
% PRE_MMN_HEP_REST_01_06_22_20220106_013811mff5ECG
% PRE_MMN_HEP_REST_01_06_22_20220106_013811mff6ECG
% PRE_MMN_HEP_REST_01_06_22_20220106_013811mff7ECG
% PRE_MMN_HEP_REST_01_06_22_20220106_013811mff8ECG
% PRE_MMN_HEP_REST_01_06_22_20220106_013811mff9ECG
% PRE_MMN_HEP_REST_01_06_22_20220106_013811mff10ECG
% PRE_MMN_HEP_REST_01_06_22_20220106_013811mff11ECG};

Total_Trials = length( ...
   all_data_ecg);
Trial_Lengths = []; %size of each trial for timekeeping
R_peaks_trial_start = [1];%keep track of trials from list of peaks
End_Trial = []; %ending indices of each trial in R_peaks_clean
Trial_Label = []; %

%%EDIT
Subject_No = 0020520220126032514; %01062220220106013811;
Condition_Matrix_Start = [2 3 4 5]; %C1 Starts @ Mat 2, C2 Starts 7 after,... until last Mat 29
%%EDIT

%% Histogram Generation %%
%% Find R_R of all peaks

%Throw-away trial for NetStation Timing
Trial_Lengths = [length(cell2mat(all_data_ecg(1)))]; %1st matrix place holder

% PER CONDITION "A"
Num_Conditions = length(Condition_Matrix_Start); % 4 conditions recorded
R_peaks = [];
for j = 1:Num_Conditions
    first = Condition_Matrix_Start(j); %matrix number per condition
    if (j==Num_Conditions) 
        last = Total_Trials; %last matrix belonging to condition
    else
        last = Condition_Matrix_Start(j+1)-1;
    end

    % PER TRIAL IN Condition "A"
    for i = first:last %for each trial
  
        ECGmat = cell2mat(all_data_ecg(i)); 

        %% filter section
        %ECGmat = Envelope for parameters?
        %OPTIONS
        %1:restructure script to analyze signal all at once and look for
        %windows of artifacts in a matrix which get converted to 0's;
        %returns the number of artifacts
        %2:use rpeaks but set parameters based on the envelope and peaks
        
        %%
        
        time_s = make_time_s(ECGmat,EEGSamplingRate); %generate time scale
        R_peaks_trial = find_Rpeaks(time_s,ECGmat); %<<---indices per trial
        
        R_peaks = [R_peaks R_peaks_trial]; %<<---indices all 

        if ~(i==Total_Trials&&j==4) %(last mat) 
            R_peaks_trial_start = [R_peaks_trial_start, length(R_peaks)+1];
        end
        
        Trial_Lengths = [Trial_Lengths length(ECGmat)]; %s.t. [applytogroup1(length mat 1), applytogroup2+gr1...]
        
        x = "error";
        if j == 1 %per condition, for each trial
             x = "MMN_";
        elseif j == 2    
             x = "ReEO";
        elseif j == 3
             x = "ReEC";
        else
             x = "HBEO";
        end
        z = repmat(x,1,length(R_peaks_trial));
        Trial_Label = [Trial_Label, z];

    end
end
%% Format Peak Via NetStation
R_peaks_NS = R_peaks;
%location in R_peaks

% Change R_peaks to include the summation of indices before it 
%NS_trial_loc = 1;
NS_sum_el = Trial_Lengths(1); %starting at mat1; %Sum of elements added to the first "n" peaks for NS

%Repeat for the number of R peaks matrices
%R_peaks_trial_start = [1<-startmat2 startmat3 startmat4....] within
%R_peaks
for i = 1:length(R_peaks_trial_start)  %there are n matrices that need to have an added length
    
    if i~=1
        NS_sum_el = Trial_Lengths(i);
    end
    
    R_peaks_NS(R_peaks_trial_start(i):end) = R_peaks_NS(R_peaks_trial_start(i):end) + NS_sum_el; %fix idx value
    R_peaks_trial_start(i);

end

%% ADDENDUM: Time for NS
R_peaks_NS; %locations of R peaks;
R_peaks_NS = (R_peaks_NS - 1); %NS adjusted peaks (for visual purposes);
time_NS = seconds(R_peaks_NS/1000); %into sec per Sample Rate = 1000Hz
time_NS_formatted = duration(time_NS,'Format','hh:mm:ss.SSS'); %formatted NS times
%% Test plot
all_ECGmat = [];
for i = 1:length(all_data_ecg)
    all_ECGmat = [all_ECGmat cell2mat(all_data_ecg(i))];
end
figure
hold on
x = 1:length(all_ECGmat);
plot(all_ECGmat);
x = R_peaks_NS;
y = ones(1,length(x)).*1000;
plot(x,y,'ro','Color','r');
hold off
%% Print NetStation Times
filename = sprintf("NS_String_R_Peaks_Subject%d_All.txt",Subject_No); 

fileID = fopen(filename,'w');
fprintf(fileID,'_%s\n',time_NS_formatted);
fclose(fileID);

%% Print all peaks
filename = sprintf("NS_R_Peaks_Subject%d_All.txt",Subject_No); 

fileID = fopen(filename,'w');
fprintf(fileID,'%d\n',R_peaks_NS);
fclose(fileID);

%% Print Labels
filename = sprintf("NS_Label_Subject%d_All.txt",Subject_No); 

fileID = fopen(filename,'w');
fprintf(fileID,'%s\n',Trial_Label);
fclose(fileID);

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
function R_index = find_Rpeaks(tm,ecgsig)
    
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
    
    y = abs(y).^2;
    %% Debug
    figure
    hold on
    plot(y)
    hold off
    %pause%  
    mean(y)

%% Filter Parameter HERE
%for post 
%MeanEnvelope; %512 resolution
mph = 1.85E5;
mpd = 0.150;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TRY TO FIND VALUES BASED ON A MEAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ENVELOPE%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


    [qrspeaks,locs] = findpeaks(y,tm,'MinPeakHeight',mph,...
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

end

%% DNC - do not change
%  Run Community Data Script to run first: /Users/samanthapitts/Desktop/Research | 
%  Mount Sinai/[Project] - HEP/Matlab File/Community Sample HEP/Figures/CommunityData12162021.m
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

close all

%From CommunityData12162021
tm = time_s;
%%%%%%%%%%
ecgsig = dataECG;
%process signal
%ecgsig = movmean(ecgsig,2);
%process signal

%%%%%%%%%%%

 

%load mit200
figure
plot(tm,ecgsig)
hold on
plot(tm,ecgsig)
xlabel('Seconds')
ylabel('Amplitude')
title('Subject - 20 - Community, Matrix 2')

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
%% Debug
figure
hold on
plot(y)
hold off
%pause%
%%

[qrspeaks,locs] = findpeaks(y,tm,'MinPeakHeight',1.75E5,...
    'MinPeakDistance',0.150); %minpeakheight 2.25E5
figure
plot(tm,y)
hold on
plot(locs,qrspeaks,'ro')
xlabel('Seconds')
title('R Peaks Localized by Wavelet Transform with Automatic Annotations')
hold off
%pause%
%% Add the expert annotations to the R-peak waveform. Automatic peak detection times are considered 
% accurate if within 150 msec of the true peak (±75 msec)
    % plot(tm(ann),y(ann),'k*')
    % title('R peaks Localized by Wavelet Transform with Expert Annotations')
%% DNC - Automated Peak Detection Compared with Raw Signal 
figure
plot(tm,ecgsig)
hold on
%edit%
locs_y = sqrt(qrspeaks); %for depiction purposes
%edit%
plot(locs,locs_y,'ro') % edit: plot(locs,qrspeaks,'ro')
xlabel('Seconds')
ylabel('Amplitude')
title('R Peaks Localized by Automatic Annotations')
hold off

%% Convert Find peaks into indices
R_index = ones(1,length(locs)).*-1;
for i = 1:length(locs)
    R_index(i) = find(tm==locs(i)); %Find index
end
%% Create a Table s.t. 
% R_index is "R Peak Indices," locs is "R Peak Times (s)," where time is calculated by 
% tm = time_s = linspace(0,EEG_ts*sample_n,sample_n); %time scale in s
varNames = {'R Peak Indices','R Peak Times (s)'}; %,'Time Vector'};
R_index = R_index';
locs = locs';
T = table(R_index,locs,'VariableNames',varNames);

varNames = {'Time Vector'};
tm = tm';
T2 = table(tm,'VariableNames',varNames);
%%edit
filename = 'R_index_Mat2_Sub20.txt';
%%edit
writetable(T,filename);


%%edit
filename = 'TimeVectorUsed_Mat2_Sub20.txt';
%%edit
writetable(T2,filename);

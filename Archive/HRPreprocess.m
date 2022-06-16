%% Visualize Data

%% FT Spectrum

%% Filter 1-40Hz

%% Visualize Filtered Signal

%% Subtract Moving Average

% Visualize Avg

% Visualize Signal without Average

%% "Zero-Phase Filtering of an Electrocardiogram Waveform" Ex
close all
%F = doFilter(ecgsig); %shift is present

figure
plot(time_s,F);
hold on
plot(time_s,ecgsig)



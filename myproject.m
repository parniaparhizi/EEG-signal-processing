% Load Cleaned EEG Dataset 
EEG = pop_loadset('filename', 'AD_preprocessed.set', ...
    'filepath', 'C:\Users\ads\Desktop\university\signal\project\SS-Project-phase1-402101396\AD\');
EEG = eeg_checkset(EEG);

% STFT Parameters
theta_band = [4 8];
gamma_band = [30 50];
fs = EEG.srate;  % Sampling rate
window_sec = 1.0;
window = hamming(round(fs * window_sec));
noverlap = round(0.9 * length(window));

% Get Trials by Event Tags (5 and 6) 
tag5_indices = find([EEG.epoch.eventtype] == 5);
tag6_indices = find([EEG.epoch.eventtype] == 6);

% Choose One Representative Channel (e.g., Cz)
chan_idx = find(strcmpi({EEG.chanlocs.labels}, 'Cz'));  % adjust if needed

% Preallocate Power Matrices
% Get time-frequency dimensions from a test trial
test_trial = squeeze(EEG.data(chan_idx,:,tag5_indices(1)));
[~, F, T, ~] = spectrogram(test_trial, window, noverlap, [], fs);
theta_idx = F >= theta_band(1) & F <= theta_band(2);
gamma_idx = F >= gamma_band(1) & F <= gamma_band(2);

power_theta_5 = zeros(length(tag5_indices), length(T));
power_gamma_5 = zeros(length(tag5_indices), length(T));
power_theta_6 = zeros(length(tag6_indices), length(T));
power_gamma_6 = zeros(length(tag6_indices), length(T));

%Compute Power Per Trial (Tag 5)
for i = 1:length(tag5_indices)
    trial = squeeze(EEG.data(chan_idx,:,tag5_indices(i)));
    [~, ~, ~, P] = spectrogram(trial, window, noverlap, [], fs);
    power_theta_5(i,:) = mean(abs(P(theta_idx,:)).^2, 1);
    power_gamma_5(i,:) = mean(abs(P(gamma_idx,:)).^2, 1);
end

% Compute Power Per Trial (Tag 6) 
for i = 1:length(tag6_indices)
    trial = squeeze(EEG.data(chan_idx,:,tag6_indices(i)));
    [~, ~, ~, P] = spectrogram(trial, window, noverlap, [], fs);
    power_theta_6(i,:) = mean(abs(P(theta_idx,:)).^2, 1);
    power_gamma_6(i,:) = mean(abs(P(gamma_idx,:)).^2, 1);
end

% Average Across Trials 
mean_theta_5 = mean(power_theta_5, 1);
mean_gamma_5 = mean(power_gamma_5, 1);
mean_theta_6 = mean(power_theta_6, 1);
mean_gamma_6 = mean(power_gamma_6, 1);

%Plot Power vs. Time Curves
figure;

subplot(2,1,1);
plot(T, mean_theta_5, 'b', 'LineWidth', 1.5); hold on;
plot(T, mean_theta_6, 'r', 'LineWidth', 1.5);
title('AD Theta Band Power (4–8 Hz)');
xlabel('Time (s)'); ylabel('Power');
legend('Chocolate', 'Rose'); grid on;

subplot(2,1,2);
plot(T, mean_gamma_5, 'b', 'LineWidth', 1.5); hold on;
plot(T, mean_gamma_6, 'r', 'LineWidth', 1.5);
title('AD Gamma Band Power (30–50 Hz)');
xlabel('Time (s)'); ylabel('Power');
legend('Chocolate', 'Rose'); grid on;

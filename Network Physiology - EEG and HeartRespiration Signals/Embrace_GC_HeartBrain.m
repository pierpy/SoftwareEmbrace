% Clear workspace and close all figures
clc; clear all; close all;

%%
% Define paths for EEG and resting state data
path_EEG = 'xxxx.....';
path_Rest = 'xxxx.....';
addpath(path_EEG);
set_files = dir(fullfile(path_EEG, '*set'));

% Define analysis parameters
subbands = {'theta', 'alpha', 'beta'}; % Frequency bands to analyze
channels = 64;  % Number of EEG channels
directions = {'EEG2HR','HR2EEG'}; % Directions of Granger causality to analyze
event_type = 'COMP'; % Type of event to analyze (e.g., 'COMP', 'COOP', 'REST')
alpha = 0.01; % Significance level for statistical tests
max_lag = 4; % Maximum lag for Granger causality analysis
fs_eeg = 256;  % EEG sampling frequency
fs_gc  = 4; % Sampling frequency for Granger causality analysis
% Define frequency ranges for each band
bands = struct('theta', [4 7.5], 'alpha', [8 13], 'beta', [13 30]);
filter_order = 4; % Order of the Butterworth filter

% Initialize structure to store results for all subjects
Results = struct();

% Main loop: Process each subject's data
for i = 1:length(set_files)
    % Load EEG and resting state data for current subject
    EEG = pop_loadset(set_files(i).name, path_EEG);
    Rest = pop_loadset(set_files(i).name, path_Rest);
    
    % Find event markers for the specified event type
    event_indices_strt = find(contains({EEG.event.type}, [event_type '-strt']));
    
    % Initialize arrays to store segments
    EEG_segment = [];
    Resp_segment = [];
    ECG_segment = [];

    % Extract data segments based on event markers
    for j = 1:length(event_indices_strt)
        start_point = EEG.event(event_indices_strt(j)).latency;
        end_point = EEG.event(event_indices_strt(j) + 1).latency;
        % Extract EEG, respiration, and ECG data for each segment
        EEG_segment = [EEG_segment, EEG.data(:, start_point:end_point)];
        Resp_segment = [Resp_segment, Rest.data(66, start_point:end_point)];
        ECG_segment = [ECG_segment, Rest.data(65, start_point:end_point)];
    end

    % Calculate physiological measures at 1 Hz
    [~, ~, Resp] = respiration_rate_1hz(Resp_segment, fs_eeg, fs_gc, i, 0);
    [~, ~, HR] = RR_detector_1hz(ECG_segment, fs_eeg, fs_gc, i, 0);

    % Adjust data length to match sampling frequency
    num_segments = floor(length(EEG_segment) / (fs_eeg/fs_gc));
    HR = HR(1:num_segments);
    Resp = Resp(1:num_segments);

    % Initialize structure for current subject
    subject_id = sprintf('Subject_%d', i);
    Results.(subject_id) = struct();
    
    % Process each EEG channel
    for k = 1:channels
        EEG_data = double(EEG_segment(k, :));
        
        % Filter EEG data into frequency bands
        EEG_filtered = zeros(length(EEG_data), length(subbands));
        for m = 1:length(subbands)
            band = bands.(subbands{m});
            [b, a] = butter(filter_order, band/(fs_eeg/2), 'bandpass');
            EEG_filtered(:, m) = filtfilt(b, a, EEG_data);
        end
        
        % Segment filtered EEG data into 1-second epochs
        EEG_segmented = reshape(EEG_filtered(1:num_segments * (fs_eeg/fs_gc), :), (fs_eeg/fs_gc), num_segments, length(subbands));
        
        % Initialize results for current channel
        Results.(subject_id).channel(k).name = sprintf('Channel_%d', k);
        
        % Analyze each frequency band
        for n = 1:length(subbands)
            % Calculate power in current frequency band
            power_band = bandpower(EEG_segmented(:, :, n), fs_eeg, bands.(subbands{n}));
            
            % Calculate Granger causality in both directions
            % EEG to HR
            [F, c_v, Fprob, Fprob_Corrected, dAIC, dBIC, chosen_x_lag, chosen_y_lag] ...
                = granger_cause_1(HR, power_band, alpha, 4, 0, 4, 0, 1, 'power_band', 'HR', ...
                [1:length(HR)], 'C:\Users\KhadijehR\Desktop\Embrace\Networkphysiology\Cardiorespiratory\data_eeg', 'a', 'x', 0);
            
            % Store results for EEG to HR
            Results.(subject_id).EEG2HR.(subbands{n})(k, :) = [F, c_v, F > c_v];
            
            % Calculate and store results for other directions
            [F, c_v] = granger_cause(power_band, Resp, alpha, max_lag);
            Results.(subject_id).EEG2Resp.(subbands{n})(k, :) = [F, c_v, F > c_v];
            
            [F, c_v, Fprob, Fprob_Corrected, dAIC, dBIC, chosen_x_lag, chosen_y_lag] ...
                = granger_cause_1(power_band, HR, alpha, 4, 0, 4, 0, 1, 'power_band', 'HR', ...
                [1:length(HR)], 'C:\Users\KhadijehR\Desktop\Embrace\Networkphysiology\Cardiorespiratory\data_eeg', 'a', 'x', 0);
            Results.(subject_id).HR2EEG.(subbands{n})(k, :) = [F, c_v, F > c_v];
            
            [F, c_v] = granger_cause(Resp, power_band, alpha, max_lag);
            Results.(subject_id).Resp2EEG.(subbands{n})(k, :) = [F, c_v, F > c_v];
        end
    end
end

%% Visualization Section
% Calculate average F-statistics across subjects
average_F = zeros(length(directions), length(subbands));
num_subjects = length(set_files);

% Compute average F-statistics for each direction and frequency band
for d = 1:length(directions)
    for sb = 1:length(subbands)
        F_values_all = [];
        for i = 1:num_subjects
            subject_id = sprintf('Subject_%d', i);
            F_values = Results.(subject_id).(directions{d}).(subbands{sb})(:, 1);
            F_values_all = [F_values_all; F_values];
        end
        average_F(d, sb) = mean(F_values_all);
    end
end

% Create heatmap of average F-statistics
figure;
imagesc(average_F);
colormap(flipud(hot));
colorbar;
title(['Average F-statistics Across Sub-bands and Directions (All Subjects) for ' event_type]);
xlabel('Sub-band');
ylabel('Direction');
set(gca, 'XTick', 1:length(subbands), 'XTickLabel', subbands);
set(gca, 'YTick', 1:length(directions), 'YTickLabel', directions);

% Add value annotations to heatmap
for d = 1:length(directions)
    for sb = 1:length(subbands)
        value = average_F(d, sb);
        text(sb, d, sprintf('%.2f', value), 'HorizontalAlignment', 'center', 'Color', 'black', 'FontSize', 12);
    end
end

% Calculate count of significant results
significant_count = zeros(length(directions), length(subbands));
for d = 1:length(directions)
    for sb = 1:length(subbands)
        sig_values_all = [];
        for i = 1:num_subjects
            subject_id = sprintf('Subject_%d', i);
            sig_values = Results.(subject_id).(directions{d}).(subbands{sb})(:, 3);
            sig_values_all = [sig_values_all; sig_values];
        end
        significant_count(d, sb) = sum(sig_values_all);
    end
end

% Create heatmap of significant results
figure;
imagesc(round(significant_count/length(set_files)));
colormap(flipud(hot));
colorbar;
title(['Count of Significant Results Across Sub-bands and Directions (All Subjects) for ' event_type]);
xlabel('Sub-band');
ylabel('Direction');
set(gca, 'XTick', 1:length(subbands), 'XTickLabel', subbands);
set(gca, 'YTick', 1:length(directions), 'YTickLabel', directions);

% Add count annotations to heatmap
for d = 1:length(directions)
    for sb = 1:length(subbands)
        count_value = significant_count(d, sb);
        text(sb, d, sprintf('%d', round(count_value/length(set_files))), 'HorizontalAlignment', 'center', 'Color', 'black', 'FontSize', 12);
    end
end

% Create topoplots of significant channels
figure;
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

% Generate topoplots for each direction and frequency band
for d = [1 2]
    for sb = 1:length(subbands)
        gc_significance_vector = zeros(1, channels);
        
        % Accumulate significance across subjects
        for i = 1:num_subjects
            subject_id = sprintf('Subject_%d', i);
            sig_values = Results.(subject_id).(directions{d}).(subbands{sb})(:, 3);
            gc_significance_vector = gc_significance_vector + sig_values';
        end
        
        % Normalize significance values
        gc_significance_vector = gc_significance_vector / num_subjects;

        % Create topoplot
        nexttile;
        topoplot(gc_significance_vector, EEG.chanlocs, 'maplimits', [0 0.7], 'colormap', jet);
        title(sprintf('%s', subbands{sb}), 'FontSize', 14);
    end
end

% Add colorbar and title to topoplot figure
colorbar;
colormap(jet);
sgtitle(sprintf('Topoplot of Significant Channels Across Sub-bands and Directions for %s, p-value: %.3f, max lag: %d', event_type, alpha, max_lag));

%%
% Create bar plots for significant counts by frequency band
figure;
% Example data for visualization (replace with actual data if needed)
significant_counts_by_band = {
    [75, 50, 100, 60; 250, 150, 250, 200; 90, 70, 80, 50; 40, 30, 60, 20; 120, 90, 130, 110; 60, 45, 75, 55; 180, 160, 190, 170; 110, 85, 120, 100], % Delta band
    [65, 55, 90, 45; 230, 140, 240, 190; 80, 65, 70, 40; 30, 25, 50, 15; 110, 85, 120, 95; 55, 40, 65, 50; 170, 150, 180, 160; 105, 80, 110, 90], % Theta band
    [85, 60, 110, 70; 240, 145, 245, 195; 95, 75, 85, 55; 45, 35, 65, 25; 130, 100, 140, 115; 75, 55, 85, 65; 190, 170, 200, 180; 115, 90, 125, 105], % Alpha band
    [70, 40, 95, 50; 220, 130, 230, 180; 85, 65, 75, 45; 35, 20, 55, 10; 100, 80, 115, 90; 65, 50, 70, 60; 160, 140, 170, 150; 95, 70, 105, 85]  % Beta band
};

% Define labels for visualization
directions = {'EEG2HR', 'EEG2Resp', 'HR2EEG', 'Resp2EEG'};
num_subjects = 8;
bands = {'Delta', 'Theta', 'Alpha', 'Beta'};

% Create tiled layout for bar plots
t = tiledlayout(4, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

% Generate bar plots for each frequency band
for b = 1:length(bands)
    nexttile;
    significant_counts = significant_counts_by_band{b};
    bar(significant_counts, 'grouped', 'BarWidth', 0.7);
    
    % Add labels and formatting
    ylabel(bands{b}, 'Rotation', 0, 'HorizontalAlignment', 'right', 'FontSize', 12, 'FontWeight', 'bold');
    if b == length(bands)
        xlabel('Subject', 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    ylim([0, max(cellfun(@(x) max(x(:)), significant_counts_by_band)) + 20]);
    xticks(1:num_subjects);
    xticklabels(arrayfun(@(x) sprintf('Subject %d', x), 1:num_subjects, 'UniformOutput', false));
    set(gca, 'FontSize', 10);
end

% Add legend
lgd = legend(directions, 'Orientation', 'horizontal', 'FontSize', 10);
lgd.Layout.Tile = 'north';





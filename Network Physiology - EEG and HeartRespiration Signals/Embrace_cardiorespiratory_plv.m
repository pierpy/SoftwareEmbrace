%% clear
clear all; close all; clc;
%%
path = 'xxx.....';
addpath(path);
mat_files     = dir(fullfile(path, '*mat'));
fs = 1024;
fs_ecg_new = 1024;    % To be fixed according to real Fs


for i=1:length(mat_files)

    load(mat_files(i).name);
    ecg_resampled = resample(double(ecg), 1,1);
    t_new = (0:length(ecg_resampled)-1)/fs_ecg_new;
  
    resp = detrend(double(resp));
    [b, a] = butter(4, 1/(fs/2), 'low');  
    filtered_resp = filtfilt(b, a, resp);
    
%     [RR_intervals, r_rest, filtered_ecg] = RR_detector2(ecg_resampled,fs_ecg_new);% if you wish to plot ecg,filtered and r peaks in a
%     %subject
    
    [resp_rest,ecg_rest,resp_comp,ecg_comp,resp_coop,ecg_coop] = segmentation(ecg_resampled,filtered_resp,trigger_info,i);
   
    results(i).plv_comp = NaN;
    results(i).plv_rest = NaN;
    results(i).plv_coop = NaN;
    results(i).hrv_coop = NaN;
    results(i).hrv_comp = NaN;
    results(i).hrv_rest = NaN;
    
    if ~isempty(ecg_rest)
        [RR_intervals, r_rest, filtered_ecg] = RR_detector2(ecg_rest, fs_ecg_new);
        results(i).plv_rest = compute_avg_plv(resp_rest, ecg_rest, r_rest, 30, fs);
        results(i).hrv_rest = mean(60./RR_intervals);
    end
    
    if ~isempty(ecg_coop)
        [RR_intervals, r_coop, filtered_ecg] = RR_detector2(ecg_coop, fs_ecg_new);
        results(i).plv_coop = compute_avg_plv(resp_coop, ecg_coop, r_coop, 30, fs);
        results(i).hrv_coop = mean(60./RR_intervals);
    end
    
    if ~isempty(ecg_comp)
        [RR_intervals, r_comp, filtered_ecg] = RR_detector2(ecg_comp, fs_ecg_new);
        results(i).plv_comp = compute_avg_plv(resp_comp, ecg_comp, r_comp, 30, fs); 
        results(i).hrv_comp = mean(60./RR_intervals);
    end

    
% Display the results
disp(['Mean PLV Comp: ', num2str(mean([results.plv_comp], 'omitnan'))]);
disp(['Mean PLV Rest: ', num2str(mean([results.plv_rest], 'omitnan'))]);
disp(['Mean PLV Coop: ', num2str(mean([results.plv_coop], 'omitnan'))]);
disp(['Mean HRV Comp: ', num2str(mean([results.hrv_comp], 'omitnan'))]);
disp(['Mean HRV Rest: ', num2str(mean([results.hrv_rest], 'omitnan'))]);
disp(['Mean HRV Coop: ', num2str(mean([results.hrv_coop], 'omitnan'))]);

%% STATISTICAL TEST
% Extract the PLV data
plv_all = [[results.plv_rest]; [results.plv_coop]; [results.plv_comp]]';
hrv_all = [[results.hrv_rest]; [results.hrv_coop]; [results.hrv_comp]]';
group = [repmat({'rest'}, 1, length([results.plv_rest])), ...
         repmat({'coop'}, 1, length([results.plv_coop])), ...
         repmat({'comp'}, 1, length([results.plv_comp]))];

% Impute missing data by replacing NaNs with the column mean
plv_all_imputed = plv_all;
hrv_all_imputed = hrv_all;
for i = 1:size(plv_all, 2)
    col_mean = nanmean(plv_all(:,i));
    plv_all_imputed(isnan(plv_all(:,i)), i) = col_mean;  % Replace NaNs with the mean
    col_mean = nanmean(hrv_all(:,i));
    hrv_all_imputed(isnan(hrv_all(:,i)), i) = col_mean;
end

[p_plv, tbl, stats] = friedman(plv_all_imputed, 1);  % Perform Friedman's test
disp(['p-value PLV: ', num2str(p_plv)]);
[p_hrv, tbl, stats] = friedman(hrv_all_imputed, 1);  % Perform Friedman's test
disp(['p-value HRV: ', num2str(p_hrv)]);

% Perform post hoc Wilcoxon signed-rank test for each pair of conditions
[p_plv_coop_rest, h] = signrank(plv_all(:,1), plv_all(:,2));  % Compare comp vs rest
[p_plv_comp_coop, h] = signrank(plv_all(:,2), plv_all(:,3));  % Compare rest vs coop
[p_plv_comp_rest, h] = signrank(plv_all(:,1), plv_all(:,3));  % Compare comp vs coop

[p_hrv_coop_rest, h] = signrank(hrv_all(:,1), hrv_all(:,2));  
[p_hrv_comp_coop, h] = signrank(hrv_all(:,2), hrv_all(:,3));  
[p_hrv_comp_rest, h] = signrank(hrv_all(:,1), hrv_all(:,3));  
%% BOXPLOT PLV
box_colors = [0.2 0.6 1; 0.8 0.4 0.1; 0.3 0.8 0.2]; % Blue, Orange, Green
figure('Units', 'inches', 'Position', [1 1 8 6]);
h = boxplot(plv_all, group, 'Widths', 0.4, 'Colors', 'k', 'Symbol', 'o');
hBox = findobj(gca, 'Tag', 'Box'); % Find all the box objects in the plot
for j = 1:length(hBox)
    patch(get(hBox(j), 'XData'), get(hBox(j), 'YData'), box_colors(mod(j-1,3)+1,:), 'FaceAlpha', 0.5); % Apply the correct color from box_colors
end

ylim([0.03, 0.41]);
yticks(0:0.05:0.4); 

% Customize plot appearance
set(gca, 'FontSize', 15, 'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out'); % Adjust font size and axis thickness
ylabel('PLV', 'FontSize', 15, 'FontWeight', 'bold');
xlabel('Condition', 'FontSize', 15, 'FontWeight', 'bold');
title('PLV Distribution Across Conditions', 'FontSize', 16, 'FontWeight', 'bold');

% Add significance marker
hold on;
y_line = 0.32; % Base height for lines
y_star_offset = 0.001; % Height offset for the stars
pvals = [p_plv_coop_rest, p_plv_comp_coop, p_plv_comp_rest];
x_pairs = [1, 2; 2, 3; 1, 3]; % X-axis locations for pairs
y_positions = [y_line, y_line + 0.02, y_line + 0.04]; % Y positions for the lines

for i = 1:length(pvals)
    if pvals(i) < 0.001
        significance = '***';  % p < 0.001
    elseif pvals(i) < 0.01
        significance = '**';   % p < 0.01
    elseif pvals(i) < 0.05
        significance = '*';    % p < 0.05
    else
        significance = 'ns';   % not significant
    end

    if ~strcmp(significance, 'ns')
        plot(x_pairs(i,:), [y_positions(i), y_positions(i)], 'k-', 'LineWidth', 1.5); % Horizontal line
        text(mean(x_pairs(i,:)), y_positions(i) + y_star_offset, significance, ...
            'FontSize', 20, 'HorizontalAlignment', 'center'); % Add significance star
    end
end

% Manually plot the individual data points using scatter
jitterAmount = 0.1; 
for i = 1:length(unique(group))
    scatter(repmat(i, size(plv_all{group == i})) + (rand(size(plv_all{group == i})) - 0.5) * jitterAmount, ...
        plv_all{group == i}, 50, box_colors(i,:), 'filled', 'MarkerFaceAlpha', 0.6); % Scatter plot for each group
end
hold off;

% Optional: Save figure in high resolution
print('PLV_Boxplot_with_Significance', '-dpng', '-r300'); % Save figure

%% BOXPLOT HRV
box_colors = [0.2 0.6 1; 0.8 0.4 0.1; 0.3 0.8 0.2]; % Blue, Orange, Green
figure('Units', 'inches', 'Position', [1 1 8 6]);
h = boxplot(hrv_all, group, 'Widths', 0.4, 'Colors', 'k', 'Symbol', 'o');
hBox = findobj(gca, 'Tag', 'Box'); % Find all the box objects in the plot
for j = 1:length(hBox)
    patch(get(hBox(j), 'XData'), get(hBox(j), 'YData'), box_colors(mod(j-1,3)+1,:), 'FaceAlpha', 0.5); % Apply the correct color from box_colors
end

ylim([50, 170]);
yticks(50:10:170); 

% Customize plot appearance
set(gca, 'FontSize', 15, 'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out'); % Adjust font size and axis thickness
ylabel('Heart Rate Variability (beats/min)', 'FontSize', 15, 'FontWeight', 'bold');
xlabel('Condition', 'FontSize', 15, 'FontWeight', 'bold');
title('Mean HRV Distribution Across Conditions', 'FontSize', 16, 'FontWeight', 'bold');

% Add significance marker
hold on;
y_line = 140; % Base height for lines
y_star_offset = 0.5; % Height offset for the stars
pvals = [p_hrv_coop_rest, p_hrv_comp_coop, p_hrv_comp_rest];
x_pairs = [1, 2; 2, 3; 1, 3]; % X-axis locations for pairs
y_positions = [y_line, y_line + 8, y_line + 16]; % Y positions for the lines

for i = 1:length(pvals)
    if pvals(i) < 0.001
        significance = '***';  % p < 0.001
    elseif pvals(i) < 0.01
        significance = '**';   % p < 0.01
    elseif pvals(i) < 0.05
        significance = '*';    % p < 0.05
    else
        significance = 'ns';   % not significant
    end

    if ~strcmp(significance, 'ns')
        plot(x_pairs(i,:), [y_positions(i), y_positions(i)], 'k-', 'LineWidth', 1.5); % Horizontal line
        text(mean(x_pairs(i,:)), y_positions(i) + y_star_offset, significance, ...
            'FontSize', 20, 'HorizontalAlignment', 'center'); % Add significance star
    end
end

% Manually plot the individual data points using scatter
jitterAmount = 0.1; 
for i = 1:length(unique(group))
    scatter(repmat(i, size(hrv_all{group == i})) + (rand(size(hrv_all{group == i})) - 0.5) * jitterAmount, ...
        hrv_all{group == i}, 50, box_colors(i,:), 'filled', 'MarkerFaceAlpha', 0.6); % Scatter plot for each group
end
hold off;



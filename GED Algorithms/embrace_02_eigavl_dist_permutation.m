%% initialization
clear; close all; clc;

path     = 'D:\PINGPONG\DATA\12_BCH1_GST2_020425\setRest\Subj2'; % replace with the path containing rest files
new_path = 'D:\PINGPONG\restSetFiles';

window_length = 256;   % in sample
num_iter      = 10000; % number of iteration for the permuation test

%% .SET files
set_files     = dir(fullfile(path, '*.set'));

%% loop over set files 
eigen_values = single(zeros(num_iter,numel(set_files)));

for data_ind = 1:numel(set_files)

    clc;

    fprintf(['WORKING ON EEG SET ' num2str(data_ind) ' / ' num2str(numel(set_files)) ' ...\n'])

    % select suitable dataset (old rest or new rest)
    data_name   = set_files(data_ind).name;
    EEG         = pop_loadset(data_name, [path '\']);

    % finding max eigenvalues
    eeg      = EEG.data;
    eig_vals = embrace_func_rest_distribution(eeg, window_length, num_iter);
    eigen_values(:,data_ind) = single(eig_vals); 
    clear EEG

    disp('|=================================================|');
end

save('eig_val_dist_clean_M','eigen_values')
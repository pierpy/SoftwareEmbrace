%% initialization
clear; close all; clc;

path = 'D:\PINGPONG\DATA\13_TTO1_SCT2_070425\setRest\Subj1'; % replace with the path containing rest files
ch = 64;
window_length   = 256;  % in sample
overlap_percent = 0;    % in percent

%% .SET files
set_files = dir(fullfile(path, '*.set'));
covR      = zeros(numel(set_files),ch,ch); % 64 is the number of EEG channels

%% loop over set files 
% numel(set_files)
for data_ind = 1:numel(set_files)

    fprintf(['WORKING ON EEG SET ' num2str(data_ind) ' / ' num2str(numel(set_files)) ' ...\n'])

    % load rest data 
    data_name    = set_files(data_ind).name;
    EEG          = pop_loadset(data_name, [path '\']);
    EEG.setname  = 'rest';
    
    % calculating covR of clean data
    eeg           = double(EEG.data);
    [reg_covR, ~] = embrace_func_compute_covR(eeg, window_length, overlap_percent);

    covR(data_ind,:,:) = reg_covR;
    
    clear ALLEEG EEG reg_covR

    fprintf('------------------------------------- \n')
end

save('covR',"covR")
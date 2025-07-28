function [covR, sample_covR] = embrace_func_compute_covR(eeg, win_length, overlap_percent)
    % eeg_data: CxT matrix containing EEG data
    % L: length of each segment (in sample)
    % overlap_percent: overlap percentage between segments
    gamma = 0.01;
    
    [C, T] = size(eeg); % Number of EEG channels and timepoints
    
    % Compute the overlap length
    overlap_length = floor(win_length * overlap_percent / 100);
    
    % Compute the number of segments
    num_segments = floor((T - overlap_length) / (win_length - overlap_length));
    
    % Compute indices for segment boundaries
    segment_starts = 1:(win_length - overlap_length):(T - win_length + 1);
    
    % Reshape EEG data into segments
    eeg_segments = zeros(C, win_length, num_segments);
    for i = 1:num_segments
        segment_start = segment_starts(i);
        segment_end = segment_start + win_length - 1;
        eeg_segments(:, :, i) = eeg(:, segment_start:segment_end);
    end
    
    % Compute covariance matrix for each segment
    sample_covR = zeros(C,C,num_segments);
    for i = 1:num_segments
        tmpd = squeeze(eeg_segments(:, :, i));
        tmpd = bsxfun(@minus,tmpd,mean(tmpd,2));
        tmp_cov = tmpd*tmpd'/(win_length);
        regularized_cov = tmp_cov*(1-gamma) + gamma*mean(eig(tmp_cov))*eye(length(tmp_cov));
        sample_covR(:, :, i) = regularized_cov;
    end
    covR = positive_definite_karcher_mean(sample_covR);
end

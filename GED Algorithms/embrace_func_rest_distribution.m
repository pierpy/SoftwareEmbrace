function eig_vals = embrace_func_rest_distribution(eeg, win_length, num_iter)
    % eeg_data: CxT matrix containing EEG data
    % L: length of each segment (in sample)
    % overlap_percent: overlap percentage between segments
    gamma = 0.01;
    
    [C, T] = size(eeg); % Number of EEG channels and timepoints
    
    % Compute the number of segments
    num_segments = floor(T / win_length);
    
    % Compute indices for segment boundaries
    segment_starts = 1:win_length:T - win_length + 1;
    
    % Reshape EEG data into segments
    eeg_segments = zeros(C, win_length, num_segments);
    for i = 1:num_segments
        segment_start = segment_starts(i);
        segment_end = segment_start + win_length - 1;
        eeg_segments(:, :, i) = eeg(:, segment_start:segment_end);
    end
    
    % Compute covariance matrix for each segment
    eig_vals = [];
    progress_bar_position = 0;
    
    for i = 1:num_iter
   
        seg1_ind = randi(num_segments);
        seg2_ind = randi(num_segments);
        if seg2_ind == seg1_ind && seg2_ind ~= 1
            seg2_ind = seg2_ind - 1;
        end
        if seg2_ind == seg1_ind && seg2_ind ~= num_segments
            seg2_ind = seg2_ind + 1;
        end

        tmpd1     = squeeze(eeg_segments(:,:,seg1_ind));
        tmpd1     = bsxfun(@minus,tmpd1,mean(tmpd1,2));
        tmp_cov  = tmpd1*tmpd1'/(win_length);
        regularized_cov1 = tmp_cov*(1-gamma) + gamma*mean(eig(tmp_cov))*eye(length(tmp_cov));

        tmpd2     = squeeze(eeg_segments(:,:,seg2_ind));
        tmpd2     = bsxfun(@minus,tmpd2,mean(tmpd2,2));
        tmp_cov  = tmpd2*tmpd2'/(win_length);
        regularized_cov2 = tmp_cov*(1-gamma) + gamma*mean(eig(tmp_cov))*eye(length(tmp_cov));

        % conventional eigen decomposition
        [~,evals] = eig(regularized_cov1,regularized_cov2);
        [evals,~] = sort(diag(evals),'descend');
        
        eig_vals = [eig_vals; evals(1)];
        
        clc;
        progress_bar_position = progress_bar_position + 1 / num_iter;
        disp(['|================= Searching Eig =================|']);
        progress_string='|';       
        for counter = 1:floor(progress_bar_position * 100 / 2)
            progress_string = [progress_string, '#'];
        end
        disp(progress_string);
        disp(['|================= ',num2str(floor(progress_bar_position * 100)),'% completed =================|']);
    end
    eig_vals = sort(eig_vals);
end

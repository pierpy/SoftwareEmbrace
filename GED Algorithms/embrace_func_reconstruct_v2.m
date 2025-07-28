function [reconstructed_eeg, num_components, eigvals_windows] = embrace_func_reconstruct_v2(eeg, reg_covR, window_length, overlap_percent, amp_threshold, eig_threshold)
    
    [C, T] = size(eeg); % Number of EEG channels and timepoints

    gamma = 0.01;
    
    % Compute the overlap length
    overlap_length = floor(window_length * overlap_percent / 100);
    
    % Compute the number of segments
    num_split = floor((T - overlap_length) / (window_length - overlap_length));
    
    % Initialize array to store filtered segments
    reconstructed_split = zeros(C, window_length, num_split);
    
    % Compute indices for segment boundaries
    split_starts = 1:(window_length - overlap_length):(T - window_length + 1);
    
    % Split EEG data into segments, filter each segment, and store filtered segments
    num_components  = zeros(4,num_split);
    eigvals_windows = zeros(64,num_split);

    for i = 1:num_split
        split_start = split_starts(i);
        split_end   = split_start + window_length - 1;
        split_data  = eeg(:, split_start:split_end);

        % ckeck if the current split is clean or not
        is_clean = embrace_func_is_clean(split_data, amp_threshold);

        if ~is_clean
        
            % reconstruction using generalized eigen decomposition
            tmpd     = split_data;
            tmpd     = bsxfun(@minus,tmpd,mean(tmpd,2));
            covS     = tmpd*tmpd'/window_length;
            reg_covS = covS*(1-gamma) + gamma*mean(eig(covS))*eye(length(covS));

            % generalized eigen decomposition
            [evecs,evals] = eig(reg_covS,reg_covR);
            [evals,sidx]  = sort(diag(evals),'descend');
            evecs         = evecs(:,sidx);
            
            % find artifactual components
            [removed_evals, no_eig, no_mad, no_knee, no_all] = embrace_func_find_noisy_components(evals, eig_threshold);

            % reconstruct clean signals by remaining components
            comp_ts       = evecs'*split_data; 
            comp_ts_trunc = comp_ts;
            comp_ts_trunc(removed_evals,:) = 0;
            eeg_reconstructed = evecs' \ comp_ts_trunc;

            % save evals and the number of rejected components
            num_components(:,i)  = [no_eig, no_mad, no_knee, no_all]';
            eigvals_windows(:,i) = evals;
        else
            eeg_reconstructed = split_data;
        end
        
        % Store the filtered segment
        reconstructed_split(:, :, i) = eeg_reconstructed;
    end
    
    % Reconstruct the main signal using the filtered segments
    reconstructed_eeg = zeros(C, T);
    for i = 1:num_split
        split_start = split_starts(i);
        split_end   = split_start + window_length - 1;
        reconstructed_eeg(:, split_start:split_end) = reconstructed_split(:, :, i);
    end

    % Make recording smooth at split time points
    filt_length  = 5;
    smoothed_eeg = reconstructed_eeg;
    temp_eeg     = movmean(reconstructed_eeg,filt_length,2);

    for i = 2:numel(split_starts)
        smoothed_eeg(:, split_starts(i) - fix(filt_length/2):split_starts(i) + fix(filt_length/2)) = ...
            temp_eeg(:, split_starts(i) - fix(filt_length/2):split_starts(i) + fix(filt_length/2));
    end
    reconstructed_eeg = smoothed_eeg;

end

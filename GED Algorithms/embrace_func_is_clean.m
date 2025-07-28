function clean_flag = embrace_func_is_clean(eeg, threshold)

    max_abs_apm = max(abs(eeg), [], "all");

    if max_abs_apm > threshold
        clean_flag = false;
    else
        clean_flag = true;
    end
    
end
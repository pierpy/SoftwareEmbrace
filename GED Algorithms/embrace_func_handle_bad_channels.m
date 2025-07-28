function [EEG,no_rmvd_chan] = embrace_func_handle_bad_channels(EEG_in)

% find bad channels by clean_rawdata
[EEG1, removed_channels] = clean_channels(EEG_in,0.5,5,[],0.5,50);

% interpolate bad channels
EEG = embrace_func_interpol(EEG1, EEG_in.chanlocs);

% number of removed channels
no_rmvd_chan = nnz(removed_channels);

end


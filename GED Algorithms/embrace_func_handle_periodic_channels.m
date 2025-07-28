function [EEG,no_periodic_channels] = embrace_func_handle_periodic_channels(EEG_in)

N  = EEG_in.pnts;
periodic_channels = false(size(EEG_in.data,1),1);

% Compute the FFT
for i = 1:size(EEG_in.data,1)
    sig = EEG_in.data(i,:);
    Y   = fft(sig);
    P2  = abs(Y/N);    
    P1  = P2(1:fix(N/2)+1); 
    P1(2:end-1) = 2*P1(2:end-1);
    
    % Analyze the spectrum to determine periodicity
    [~, peakFreqIdx]  = max(P1);           
    total_power       = sum(P1.^2);        
    peak_power        = P1(peakFreqIdx)^2; 
    
    % Calculate the ratio of the peak power to the total power
    power_ratio = peak_power / total_power;
    
    % Define a threshold for "pure" periodicity, e.g., 40% of power concentrated at one frequency
    if power_ratio > 0.4
        periodic_channels(i,1) = true;
    end
end

EEG          = EEG_in;
EEG.data     = EEG.data(~periodic_channels,:);
EEG.nbchan   = size(EEG.data, 1);
EEG.chanlocs = EEG.chanlocs(~periodic_channels);
EEG          = embrace_func_interpol(EEG, EEG_in.chanlocs);

no_periodic_channels = sum(periodic_channels);

end
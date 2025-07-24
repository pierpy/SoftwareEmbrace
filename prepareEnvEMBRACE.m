function envOut = prepareEnvEMBRACE(emg, frequency)

if size(emg,2)>size(emg,1)
    emg = emg';
end

emg_filt = zeros(size(emg));
for i=1:size(emg,2)
    [b,a] = butter(3, [20, 450]/frequency, 'bandpass');
    emg_filt(:,i) = filtfilt(b, a, emg(:,i));
    [b,a] = butter(3, [49.5, 50.5]/frequency, 'stop');
    emg_filt(:,i) = filtfilt(b, a, emg_filt(:,i));
    [b,a] = butter(3,5/frequency,'low');
    env_tmp(:,i) = filtfilt(b, a, abs(emg_filt(:,i)));
end
env_tmp(env_tmp<(1e-4)*max(env_tmp(:))) = (1e-4)*max(env_tmp(:));

envOut = env_tmp;

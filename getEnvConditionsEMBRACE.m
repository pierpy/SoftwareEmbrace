function envOut = getEnvConditionsEMBRACE(condition)

data_folder = uigetdir();
data_folder = fullfile(data_folder, "motion+emg");

files = dir(strcat(data_folder, "/*_", condition, "*.tdf"));
env = [];
for file=1:length(files)
    try
        data_file = fullfile(files(file).folder, files(file).name);
        [~,frequency,~,~,emgData] = tdfReadDataEmg(data_file);
        env_tmp = prepareEnvEMBRACE(emgData, frequency);
        env = [env;env_tmp];
    catch
        disp(strcat("Error in the file: ", files(file).name));
    end
end

for i = 1:size(env,2)
    env(:,2) = env(:,2)/prctile(env(:,2),95);
end

envOut = env;
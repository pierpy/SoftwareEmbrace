clc
clear
close all

% Defines the paths to the data.
config.path.eeg  = '../../data/eeglab/';
config.path.sens = '../../data/segments/';
config.path.patt = '*_*';


% Creates the output folder, if required.
if ~exist ( config.path.sens, 'dir' ), mkdir ( config.path.sens ), end


% Lists the dyads.
dyads = dir ( sprintf ( '%s%s', config.path.eeg, config.path.patt ) );

% Goes through all the dyads.
for dindex = 1: numel ( dyads )

    % Gets the dyad.
    dyad  = dyads ( dindex ).name;

    % Gets the subject code.
    hits  = regexp ( dyad, '^\d+_([A-Z]+1_[A-Z]+2)_', 'tokens' );
    subj  = hits {1} {1};
    subj  = strrep ( subj, '_', '-' );


    % Lists the files for this dyad.
    files = dir ( sprintf ( '%s%s/setTaskClean/Subj1/*.set', config.path.eeg, dyad ) );
    
    % Gets the number of files.
    nfile = numel ( files );

    % Iitializes the trial data structure.
    trialdata           = [];
    trialdata.label     = {};
    trialdata.time      = cell ( nfile, 1 );
    trialdata.trial     = cell ( nfile, 1 );
    trialdata.fsample   = [];
    trialdata.trialinfo = nan ( nfile, 3 );

    % Initializes the list of bad channels.
    bads  = {};

    % Goes through each file.
    for findex = 1: numel ( files )
        
        % Lists the data and metadata files.
        meta1 = sprintf ( '%s%s/setTaskClean/Subj1/%s', config.path.eeg, dyad, files ( findex ).name );
        data1 = regexprep ( meta1, '.set$', '.fdt' );
        meta2 = sprintf ( '%s%s/setTaskClean/Subj2/%s', config.path.eeg, dyad, files ( findex ).name );
        data2 = regexprep ( meta2, '.set$', '.fdt' );

        % Determines the condition.
        hit   = regexp ( files ( findex ).name, '(Comp|Coop|SOLO)', 'tokens' );
        cond1 = hit {1} {1};
        cond1 = find ( strcmp ( cond1, { 'Comp' 'Coop' 'SOLO' } ) );
        hit   = regexp ( files ( findex ).name, '(Easy|Hard|SUBJ1|SUBJ2)', 'tokens' );
        cond2 = hit {1} {1};
        cond2 = find ( strcmp ( cond2, { 'Easy' 'Hard' 'SUBJ1' 'SUBJ2' } ) );

        % Loads the metadata.
        meta1 = load ( '-mat', meta1 );
        meta2 = load ( '-mat', meta2 );
        
        % Gets the EEG labels.
        label1 = strtrim ( { meta1.chanlocs.labels }' );
        label1 = strcat ( 'A_', label1 );
        label2 = strtrim ( { meta2.chanlocs.labels }' );
        label2 = strcat ( 'B_', label2 );

        % Gets the list of bad channels.
        bads1  = strtrim ( { meta1.chaninfo.removedchans.labels }' );
        bads1  = strcat ( 'A_', bads1 );
        bads2  = strtrim ( { meta2.chaninfo.removedchans.labels }' );
        bads2  = strcat ( 'B_', bads2 );
        
        % Gets the time vectors.
        time1  = meta1.times * 1e-3;
        time2  = meta2.times * 1e-3;
        % nsamp  = min ( numel ( time1 ), numel ( time2 ) );
        nsamp  = 19106;
        
        % % Checks that the times are equal.
        % if ~isequal ( time1, time2 )
        %     warning ( 'Time vectors are not equal (%.i vs. %i).', numel ( time1 ), numel ( time2 ) )
        % end

        % Keeps only the desired time.
        time1  = time1 ( 1: nsamp );
        time2  = time2 ( 1: nsamp );
        
        % Gets the sampling rates.
        srate1 = meta1.srate;
        srate2 = meta2.srate;
        
        % Checks that the times are equal.
        if srate1 ~= srate2
            error ( 'Sampling rates are not equal.' )
        end
        
        
        % Loads the data.
        fid = fopen ( data1, 'rb' );
        data1 = fread ( fid, [ meta1.nbchan nsamp ], 'float' );
        fclose ( fid );
        fid = fopen ( data2, 'rb' );
        data2 = fread ( fid, [ meta2.nbchan nsamp ], 'float' );
        fclose ( fid );
        
        % Concatenates the labels and the data.
        label  = cat ( 1, label1, label2 );
        data   = cat ( 1, data1, data2 );
        bads   = union ( bads, cat ( 1, bads1, bads2 ) );
        
        % Unifies the sampling rate and time.
        time   = time1;
        srate  = srate1;
        
        
        % Stores the data and time.
        trialdata.time  { findex } = time;
        trialdata.trial { findex } = data;

        % Stores the condition.
        trialdata.trialinfo ( findex, : ) = [ findex cond1 cond2 ];
        
        % Adds the sampling rate and channel labels, if required.
        if findex == 1
            trialdata.label   = label;
            trialdata.fsample = srate;
        
        % Checks the consistence of the sampling rate and channel labels.
        else
            if ~isequal ( trialdata.label, label )
                error ( 'Inconsistent channels.' )
            end
            if trialdata.fsample ~= srate
                error ( 'Inconsistent sampling rate.' )
            end
        end
    end

    % Sets the trial information.
    trialinfo           = [];
    trialinfo.trialfile = trialdata.trialinfo ( :, 1 );
    trialinfo.trialdef  = trialdata.trialinfo;
    trialinfo.trialpad  = 2 * ones ( nfile, 1 );

    % Sets the channel information.
    chaninfo            = [];
    chaninfo.bad        = bads;

    % Compiles the information.
    sensdata            = [];
    sensdata.subject    = subj;
    sensdata.task       = 'all';
    sensdata.stage      = '';
    sensdata.channel    = 'EEG';
    sensdata.trialdata  = trialdata;
    sensdata.chaninfo   = chaninfo;
    sensdata.trialinfo  = trialinfo;

    % Saves the data.
    save ( '-v6', sprintf ( '%s%s_%s%s_%s.mat', config.path.sens, sensdata.subject, sensdata.task, sensdata.stage, sensdata.channel ), '-struct', 'sensdata' )
end

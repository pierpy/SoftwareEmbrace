function erfdata = myft_timelock ( cfg, trialdata )

% Sanitizes the input.
if ~isfield ( cfg, 'latency'    ), cfg.latency    = [ -Inf +Inf ];                     end
if ~isfield ( cfg, 'baseline'   ), cfg.baseline   = [ -Inf +0.0 ];                     end
if ~isfield ( cfg, 'channel'    ), cfg.channel    = trialdata.label;                   end
if ~isfield ( cfg, 'trials'     ), cfg.trials     = true ( size ( trialdata.trial ) ); end
if ~isfield ( cfg, 'baseline'   ), cfg.baseline   = [];                                end
if ~isfield ( cfg, 'keeptrials' ), cfg.keeptrials = false;                             end

if ischar ( cfg.latency ) && strcmpi ( cfg.latency, 'all' ), cfg.latency = [ -Inf +Inf ];                     end
if ischar ( cfg.trials  ) && strcmpi ( cfg.trials,  'all' ), cfg.trials  = true ( size ( trialdata.trial ) ); end

if ischar ( cfg.keeptrials ), cfg.keeptrials = strcmpi ( cfg.keeptrials, 'yes' ); end


% Extracts the data and the time vector.
trial   = cat ( 3, trialdata.trial { cfg.trials } );
time    = trialdata.time {1};


% Gets the indexes of the desired channels.
label   = ft_channelselection ( cfg.channel, trialdata.label );
chindex = ismember ( trialdata.label, label );

% Keeps only the desired channels.
trial   = trial ( chindex, :, : );
label   = trialdata.label ( chindex );


% Gets the indexes of the data and baseline samples.
lindex  = time >= cfg.latency  (1) & time <= cfg.latency  (2);
bindex  = time >= cfg.baseline (1) & time <= cfg.baseline (2);

% Gets the baseline.
bline   = trial ( :, bindex, : );

% Gets the data.
trial   = trial ( :, lindex, : );
time    = time ( lindex );


% Corrects the data by the baseline, if requested.
if numel ( bline )
    trial   = bsxfun ( @minus, trial, mean ( bline, 2 ) );
    bline   = bsxfun ( @minus, bline, mean ( bline, 2 ) );
end


% Calculates the average and variance of the ERF.
erfavg  = mean ( trial, 3 );
erfvar  = var  ( trial, [], 3 );


% Calculates the data and baseline covariances.
trialcov = trial ( :, : ) * trial ( :, : )' / size ( trial, 3 ) / ( size ( trial, 2 ) - 1 );
blinecov = bline ( :, : ) * bline ( :, : )' / size ( bline, 3 ) / ( size ( bline, 2 ) - 1 );


if ~isfield ( trialdata, 'elec' )
    trialdata.elec = [];
end

% Creates the FT timelock structure.
if cfg.keeptrials
    erfdata           = [];
    erfdata.label     = label;
    erfdata.trial     = permute ( trial, [ 3 1 2 ] );
    erfdata.avg       = erfavg;
    erfdata.var       = erfvar;
    erfdata.time      = time;
    erfdata.dimord    = 'rpt_chan_time';
    erfdata.dof       = size ( trial, 3 ) * ones ( size ( trial, 1 ), size ( trial, 2 ) );
    erfdata.trialinfo = trialdata.trialinfo;
    if isfield ( erfdata, 'elec' )
        erfdata.elec      = trialdata.elec;
    end
    if isfield ( erfdata, 'grad' )
        erfdata.grad      = trialdata.grad;
    end
    erfdata.cov       = trialcov;
    erfdata.noisecov  = blinecov;
else
    erfdata           = [];
    erfdata.label     = label;
    erfdata.avg       = erfavg;
    erfdata.var       = erfvar;
    erfdata.time      = time;
    erfdata.dimord    = 'chan_time';
    erfdata.dof       = size ( trial, 3 ) * ones ( size ( trial, 1 ), size ( trial, 2 ) );
    if isfield ( erfdata, 'elec' )
        erfdata.elec      = trialdata.elec;
    end
    if isfield ( erfdata, 'grad' )
        erfdata.grad      = trialdata.grad;
    end
    erfdata.cov       = trialcov;
    erfdata.noisecov  = blinecov;
end

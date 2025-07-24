clc
clear
close all

% config.path.segs = '../../data/segments/';
config.path.pow  = '../../data/spectra/dpss_05_new/';
config.path.figs = '../../figs/spectra/dpss_05_new/';
config.path.patt = '*.mat';

% Selects which versions of the figure to save.
config.savefig       = false;

% Selects the layout to plot the data.
config.layout        = 'ant64dry_hyper.mat';
config.occipital     = { '5*' '7*' '8*' '9*' '10*' };


% Adds the functions folders to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );

% Adds, if needed, the FieldTrip folder to the path.
myft_path


% Generates the output folder, if needed.
if ~exist ( config.path.figs, 'dir' ), mkdir ( config.path.figs ); end


% Generates the layout.
layout   = ft_prepare_layout ( config );

% Gets the list of files.
files    = dir ( sprintf ( '%s%s', config.path.pow, config.path.patt ) );

% Goes through each file.
for findex = 1: numel ( files )
    
    % Loads the data.
    powdata   = load ( sprintf ( '%s%s', config.path.pow, files ( findex ).name ) );
    freqdata  = powdata.freqdata;
    
    
    figure ( 'Position', [ 105 558 1135 420 ] )
    
    % Draws the spectra.
    cfg           = [];
    cfg.layout    = layout;
    cfg.comment   = ' ';
    
    subplot ( 2, 2, [ 1 3 ] )
    ft_multiplotER ( cfg, freqdata )
    
    
    % Draws the average occipital spectrum.
    cfg           = [];
    cfg.layout    = layout;
    cfg.channel   = config.occipital;
    
    subplot ( 2, 2, 2 )
    ft_singleplotER ( cfg, freqdata )
    title ''
    
    % Draws the alpha topology.
    cfg           = [];
    cfg.layout    = layout;
    cfg.xlim      = [ 7 12 ];
    cfg.comment   = 'no';
    
    subplot ( 2, 2, 4 )
    ft_topoplotER ( cfg, freqdata )
    
    fig = gcf;
    fig.PaperOrientation = 'portrait';
    
    % Saves the figure.
    print ( '-dpng', sprintf ( '%s%s_%s%s_%s.png', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    
    if config.savefig
        savefig ( sprintf ( '%s%s_%s%s_%s.png', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    end
    
    close all
end

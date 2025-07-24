#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@author: Ricardo Bruña

"""

import re
import datetime

import mne
import numpy


# Lists the valid MNE objects.
mnevalid = (
    mne.io.BaseRaw,
    mne.BaseEpochs )

# Sets the verbosity level for MNE.
mne.set_log_level ( verbose = 'ERROR' )


def build_raw ( info, data, montage = None ):
    
    
    # Lists the channels in the data.
    ch_label = info [ 'channels' ] [ 'label' ]
    
    
    # If no montage assumes the standard 10-05.
    if montage is None:
        montage  = mne.channels.make_standard_montage ( 'standard_1005' )
    
    
    # Identifies the EEG, EOG, ECG, and EMG channels.
    ind_eeg  = numpy.where ( numpy.in1d ( ch_label, montage.ch_names ) )
    ind_eog  = numpy.where ( [ re.search ( 'EOG', label ) != None for label in ch_label ] )
    ind_ecg  = numpy.where ( [ re.search ( 'CLAV', label ) != None for label in ch_label ] )
    ind_emg  = numpy.where ( [ re.search ( 'EMG', label ) != None for label in ch_label ] )
    
    # Marks all the channels as EEG.
    ch_types = numpy.array ( [ 'eeg' ] * len ( ch_label ) )

    # Sets the channel types.
    ch_types [ ind_eeg ] = 'eeg'
    ch_types [ ind_eog ] = 'eog'
    ch_types [ ind_ecg ] = 'ecg'
    ch_types [ ind_emg ] = 'emg'
    
    
    # Creates the MNE-Python information object.
    mneinfo  = mne.create_info ( 
        ch_names  = list ( info [ 'channels' ] [ 'label' ] ),
        sfreq     = info [ 'sample_rate' ],
        ch_types  = list ( ch_types ) )
    
    # Adds the montage, if provided.
    if montage is not None:
        mneinfo.set_montage ( montage )
    
    
    # Creates the MNE-Python raw data object.
    mneraw   = mne.io.RawArray ( data.T, mneinfo, verbose = False )
    
    # Overwrites the default parameters.
    mneraw.set_meas_date ( info [ 'acquisition_time' ] )
    
    # Adds the calibration factor.
    mneraw._cals = numpy.ones ( len ( ch_label ) )
    
    # Marks the 'active' channels.
    mneraw._read_picks = [ numpy.arange ( len ( ch_label ) ) ]
    
    
    # Gets the information about the impedances, if any.
    if 'impedances' in info:
        
        # Takes only the first measurement.
        if len(info['impedances']) > 0:
            impmeta    = info [ 'impedances' ] [0]
            impedances = impmeta [ 'measurement' ]

            # Fills the extra information for MNE.
            for channel, value in impedances.items ():

                impedances [ channel ] = {
                    'imp':           value,
                    'imp_unit':      impmeta [ 'unit' ],
                    'imp_meas_time': datetime.datetime.fromtimestamp ( impmeta [ 'time' ] ) }


            # Adds the impedances to the MNE object.
            mneraw.impedances = impedances
    
    
    # Gets the annotations, if any.
    annotations = mne.Annotations (
        [ annot [ 'onset' ] for annot in info [ 'events' ] ],
        [ annot [ 'duration' ] for annot in info [ 'events' ] ],
        [ annot [ 'description' ] for annot in info [ 'events' ] ] )
    
    # Adds the annotations to the MNE object.
    mneraw.set_annotations ( annotations )
        
    
    # Returns the MNE Raw object.
    return mneraw

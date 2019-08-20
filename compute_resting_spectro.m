%%compute_resting_spectro(EEG)
%
%This function computes spectral features, channel x time densities, absolute and relative power for fixed frequency bands, individual alpha peak, and individual frequency bands (based on alpha peak). 
%
% INPUT:
%  EEG: EEG structure with continuous 105-channel EEG (EGI system Langer-Lab cap layout) data.
%       The data should be based on eyes open and eyes closed from resting state recordings
%       EEG.events structure should contain the following types:
%       1 x '90' event (start of the recording)
%       5 x '20' events (eyes open condition, 20 sec duration)
%       5 x '30' events (eyes closed condition, 40 sec duration)
%       The data should already be preprocessed (e.g. using the Automagic toolbox) and 
%       should not contain eye, ecg, or other than EEG channels
%
% OUTPUT:
%   EEG: EEG structure, EEG.data contains notch (48-52 Hz) and bandpass (1-90 Hz) filtered continuous data
%        EEG.spectro contains several substructures containing the outputs:
%        eyesclosed / eyesopen : refers to the different experimental condition
%        welch / fft : refers to different methods used for spectral decomposition
%        specdata : channel x frequency matrix of power spectral densities 
%        freqs : frequency x 1 vector of frequencies
%        fbands : structure containing absolute and relative power for fixed frequency bands (delta to gamma),
%                 across all channels and for selected frequency clusters
%        alphaPeak : structure containing individual alpha peak frequency and amplitude using different detection methods 
%                    (Maximum, Derivative, Amplitude)
%        indfbands : individual frequency bands, defined relative to individual alpha peak
%        (see subfunctions for detailled description of how outputs are computed)
%
% DEPENDENCIES:
%   eeglab added to the matlab path
%   script requires the current directory to be set to the location of the script location (as it needs to add the function folder to the path)
%
%
% christian.pfeiffer@uzh.ch
% 20.08.2019
%
function EEG = compute_resting_spectro(EEG)

  %loads processing settings and adds 'function' folder to matlab path
  settings = locfun_default_settings();

  %notch filter
  EEG = pop_eegfiltnew(EEG, settings.spectro.notch.lpf,    settings.spectro.notch.hpf,   [], 1); %[]=default filter order, and 'revfilt'=1 for notch
  
  %bandpass filter
  EEG = pop_eegfiltnew(EEG, settings.spectro.bandpass.lpf, settings.spectro.bandpass.hpf      );

  %re-referencing to average
  if settings.averageref == 1
      EEG = pop_reref(EEG,[]);
  end

  %checks for the presence of the correct number of events (i.e. 5 x '20' and 5 x '30' trigger)
  %and remove leading or trailing whitespaces, and additional events of no interest
  [EEG,isvalid] = fix_event_structure(EEG);

  %if the event structure does not contain a valid amount of events, don't continue
  if ~isvalid

    disp(['..skipping file. Number of events is not valid.'])

  %otherwise, continue the analysis..
  else

    %% SPECTRO ANALYSIS 
    for eyes = {'eyesclosed','eyesopen'}
      
      %Segmentation, cuts out segments interest (e.g. eyesclosed only), and concatenates them into a "continuous" dataset
      EEGtmp = RestingSegment(EEG,settings.segment.(eyes{1}));
      
      %Spectrogram (using spectopo) for all good segments of the data (i.e. 1min30sec for eyesopen and 3min10sec for eyesclosed)
      EEG.spectro.(eyes{1}) = RestingSpectro(EEGtmp, settings.spectro, eyes{1});
      
      clear EEGtmp;
        
    end

    %% INDIVIDUAL ALPHA PEAK
    for eyes = {'eyesclosed','eyesopen'}
      for mn = {'fft','welch'}
        
        %Individual alpha peak for all the data
        EEG.spectro.(eyes{1}).(mn{1}).alphaPeak = RestingAlphaPeak(EEG.spectro.(eyes{1}).(mn{1}), settings, eyes{1});
        
      end
    end

    %% Individually defined frequency bands and frequency ratios
    for eyes = {'eyesclosed','eyesopen'}
      for mn = {'fft','welch'}
      
        %for all the data
        EEG.spectro.(eyes{1}).(mn{1}) = IndividualSpectro(EEG.spectro.(eyes{1}).(mn{1}),settings.spectro);

      end
    end

  end
end


%settings = locfun_default_settings()
%
%This function returns a structure 'settings' with parameters for resting data processing
%and adds the folder 'function' to the matlab path
%
function settings = locfun_default_settings()

  settings = struct();

  addpath([pwd,filesep,'functions',filesep])

  %% Plotting settings -----------------------------------------------
  settings.figure.visible = 'on';

  %% Eye tracking settings -------------------------------------------
  settings.ET_resting.trig_eyeo = 20;
  settings.ET_resting.trig_eyec = 30;
  settings.ET_resting.seg_time = 500; % how much should be cutted in addition on the beginning and end of each segment in ms
  settings.ET_resting.outlierstd = 2;

  %% do average re-referencing
  settings.averageref = 1; 

  %% Segmentation 
  settings.segment = {};
  settings.segment.fun = 'restingsegment';
  settings.segment.path = {};
  settings.segment.eyesclosed.events = '30'; % == eyes closed
  settings.segment.eyesclosed.timelimits = [1 39]; % cut out 1 sec at the onset and 1 sec before the end
  settings.segment.eyesopen.events = '20'; % == eyes open
  settings.segment.eyesopen.timelimits = [1 19]; % cut out 1 sec at the onset and 1 sec before the end

  %% spectrogram Analysis 
  settings.spectro = {};
  settings.spectro.fun = 'restingspectro';
  settings.spectro.path = {};
  settings.spectro.notch.lpf = 48; %
  settings.spectro.notch.hpf = 52; %
  settings.spectro.bandpass.lpf = 1; %
  settings.spectro.bandpass.hpf = 90; %
  settings.spectro.winlength = 1000; % = 2 Seconds
  settings.spectro.timelimits = [0 1000]; % 0 to 2 Seconds
  settings.spectro.mvmax = 90; % maximum millivoltage to clean data
  settings.spectro.numsegments = [5 15 30 45 60 75 90]; %how many good segements to include for separate analyses
  settings.spectro.fbands = {};
  settings.spectro.doplot= 1;

  % the frequencies of interest. Define the lower and upper limits of the
  % relative power normalization 
  fbands =    {'delta','theta','alpha','alpha1','alpha2','beta','beta1','beta2','beta3','gamma','gamma1','gamma2'};
  lowfreqs =  [  1.5,   4.0,     8.5,    8.5,    10.5,    12.5,  12.5,   18.5,   21.5,   30.5,   30.5,    45.5   ];
  highfreqs = [  3.5,   8.0,    12.0,   10.0,    12.0,    30.0,  18.0,   21.0,   30.0,   80.0,   45.0,    80.0   ];

  for i=1:length(fbands)
      settings.spectro.fbands(i).name = fbands{i};
      settings.spectro.fbands(i).lowfreqs = lowfreqs(i);
      settings.spectro.fbands(i).highfreqs = highfreqs(i);
  end
  clear fbands lowfreqs highfreqs;

  % electrode arrays to average frequency bands
  eleclusters.names = {'l_front','m_front','r_front', 'l_pari','m_pari','r_pari'};
  eleclusters.chans = {{'E33' , 'E26' , 'E22' , 'E34' , 'E27' , 'E23' , 'E35' , 'E28' , 'E24' , 'E19' , 'E36' , 'E29'  , 'E20' , 'E30' , 'E13' }, ...
                       {'E18' , 'E12' , 'E6'  , 'E7'  , 'E31' , 'E15' , 'E16' , 'E11' , 'Cz'  , 'E10' , 'E5'  , 'E106' , 'E80' }, ...
                       {'E9'  , 'E4'  , 'E118', 'E112', 'E105', 'E3'  , 'E124', 'E111', 'E104', 'E2'  , 'E123', 'E117' , 'E110', 'E116', 'E122'}, ...
                       {'E45' , 'E50' , 'E58' , 'E65' , 'E70' , 'E46' , 'E51' , 'E59' , 'E66' , 'E41' , 'E47' , 'E52'  , 'E60' , 'E42' , 'E53' , 'E37' }, ... 
                       {'E54' , 'E61' , 'E67' , 'E71' , 'E75' , 'E55' , 'E62' , 'E72' , 'E79' , 'E78' , 'E77' , 'E76'  }, ...
                       {'E83' , 'E90' , 'E96' , 'E101', 'E108', 'E84' , 'E91' , 'E97' , 'E102', 'E85' , 'E92' , 'E98'  , 'E103', 'E86' , 'E93' , 'E87' }};

  for i=1:length(eleclusters.names)
      settings.spectro.eleclusters(i).names = eleclusters.names{i};
      settings.spectro.eleclusters(i).chans = eleclusters.chans{i};
  end                           
               
  %% settings for alpha peak
  %  %  Reference: Grandy et al. 2014: Mean spectrum of posterior electrodes
  %  1. Alpha individual peak = largest power between 7.5 and 12.5
  %  2. Weighted mean: IAF = (Sum(a(f) x  f))/(Sum a(f)). (Klimesch)
  %  3. First derivative changeing point

  % Alpha amplitude was defined as the mean amplitude
  % of the frequency spectrum of the 17 posterior electrodes
  % 1Hz around the IAF. 

  settings.alphapeak.postelectrodes = {   'E53' , 'E61' , 'E62' , 'E78' , 'E86' , 'E52' , 'E60' , 'E67' , ...
                                          'E72' , 'E77' , 'E85' , 'E92' , 'E59' , 'E66' , 'E71' , 'E76' , ...
                                          'E84' , 'E91' , 'E70' , 'E75' , 'E83' };

  settings.alphapeak.type = 'deriv'; % 'max', 'wmean'
  settings.alphapeak.lower = 7.5; %% reference: Grandy et al., 2014 use 'deriv'
  settings.alphapeak.upper = 12.5;
  settings.alphapeak.window = 1; % Amplitude +- 1 Hz around peak is the mean individual alpha amplitude

end
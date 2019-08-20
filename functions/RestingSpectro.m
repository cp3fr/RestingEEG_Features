% Calculates the spectral density of the Resting EEG 
% EEG = eeglab EEG
% settings = settings from RestingCreateSettings.
% eyes = 'eyesclosed', 'eyesopen'
% output: - EEG.specdata  Averaged power spectral density (PSD)
%         - EEG.freqs = frequency bins
% chanlocs = channel locations (for plotting)
% name = method of frequency decomposition (welch, ftt)
% eyes = see above

function out = RestingSpectro(EEG,settings,eyes)
fprintf('\n:::RestingSpectro...\n')

%% extract snippets of X seconds out of EEG.data
%segments the data in non-overlapping epochs of 2sec/1000sp length
EEG.data = epoch(EEG.data,1:settings.winlength:size(EEG.data,2),settings.timelimits);
EEG.pnts = size(EEG.data,2);
EEG.times = 1: 1000/EEG.srate: size(EEG.data,2)*1000/EEG.srate;
EEG.trials = size(EEG.data,3);

%% find the good segments 
% using <90uV amplitude criterion
% this could be modified with other criteria...
EEG.goodsegments = [];
for g=1:size(EEG.data,3)
    A = EEG.data(:,:,g);
    if all(abs(A(:)) < settings.mvmax)  % tests if all recorded samples are below mV max (e.g. 90 mV)
        EEG.goodsegments(end+1) = g;
    end
end

        %% do the pwelch spectrogram using spectopo function
        % Not sure if these are the optimal settings for spectopo?
        [EEG.welch.specdata, EEG.welch.freqs] = spectopo(EEG.data(:,:,EEG.goodsegments),0,EEG.srate,'freqfac',2,'plot','off');
        EEG.welch.specdata = db2pow(EEG.welch.specdata);
        %imagesc(EEG.welch.specdata)
        %plot(mean(EEG.welch.specdata(:,1:100)))
        
        EEG.welch.fbands = computeFbands(EEG.welch.specdata,EEG.welch.freqs,settings,EEG.chanlocs);
 
        %% FFT is calculated
        EEG.fft.freqs = (0:150)*EEG.srate/settings.winlength;
        spec_p=zeros(EEG.nbchan,length(EEG.fft.freqs));
                   
        for g = EEG.goodsegments
            temp = abs(fft(squeeze(EEG.data(:,:,g)),[],2))./(settings.winlength/2); % the fft is applied to each channel in each epoch
            spec_p = spec_p + temp(:,1:length(EEG.fft.freqs));
        end
        EEG.fft.specdata = spec_p/size(EEG.data,3);


        EEG.fft.fbands = computeFbands(EEG.fft.specdata,EEG.fft.freqs,settings,EEG.chanlocs);
    

%% compute the frequency bands
    function fbands = computeFbands(specdata,freqs,settings,chanlocs)
                
        for f = 1:length(settings.fbands)
            % here we could add a check for filters!
            ind_notch = freqs >= settings.notch.lpf & freqs <= settings.notch.hpf;
            ind_band = freqs >= settings.fbands(f).lowfreqs & freqs <= settings.fbands(f).highfreqs;
            ind_band_no_notch = ind_band & ~ind_notch;
            ind_range = freqs >= settings.fbands(1).lowfreqs & freqs <= settings.fbands(end).highfreqs;
            ind_range_no_notch = ind_range & ~ind_notch;
            
            fbands(f).name = settings.fbands(f).name;
            fbands(f).lowfreqs = settings.fbands(f).lowfreqs;
            fbands(f).highfreqs = settings.fbands(f).highfreqs;
            fbands(f).absmean = nanmean(specdata(:,ind_band_no_notch),2);
            fbands(f).relmean = nanmean(specdata(:,ind_band_no_notch),2) ./ nanmean(specdata(:,ind_range_no_notch),2);
            
            % average over electrode clusters
            for  k = 1:length(settings.eleclusters)
                clusterindex = ismember({chanlocs.labels}, settings.eleclusters(k).chans);
                fbands(f).elecluster(k).names = settings.eleclusters(k).names;
                fbands(f).elecluster(k).absmean = nanmean(fbands(f).absmean(clusterindex));
                fbands(f).elecluster(k).relmean = nanmean(fbands(f).relmean(clusterindex));
            end
        end
              
    end

%% cleanup
EEG.pnts = size(EEG.data,2);
EEG.times = 1:1000/EEG.srate: size(EEG.data,2)*1000/EEG.srate;
EEG.trials = size(EEG.data,3);

%% Save output
out.welch = EEG.welch;
out.fft = EEG.fft;
out.goodsegments = EEG.goodsegments;
out.fft.chanlocs = EEG.chanlocs;
out.welch.chanlocs = EEG.chanlocs;
out.fft.name = 'fft';
out.welch.name = 'welch';
out.fft.eyes = eyes;
out.welch.eyes = eyes;

end
% Extracts the Alpha Peak 
% EEG = eeglab EEG
% settings = settings from RestingCreateSettings
% eyes = 'eyesclosed', 'eyesopen'
% output: 
% alphapeakfreqMax: find the maximimum point in the search range 
% alphapeakamplitudeMax: the amplitude at this point
% alphapeakfreqDerivative: first derivative crossing zero from left side 
% alphapeakamplitudeDerivative: the amplitude at this point
% alphapeakfreqGravity: 'gravity method' weighted mean
% alphapeakamplitudeGravity: the amplitude at this point
% uniquelyIdentifiable: [0/1] is there one peak identifiable with max and deriv.
% method? 


function out = RestingAlphaPeak(EEG,settings,eyes)
fprintf('\n:::RestingAlphaPeak...\n')

% define the window of interest
idxlower = find(EEG.freqs == settings.alphapeak.lower);
idxupper = find(EEG.freqs == settings.alphapeak.upper);
freqRes = EEG.freqs(2)-EEG.freqs(1);
freqRange =[settings.alphapeak.lower:freqRes:settings.alphapeak.upper];

%%  average the posterior cluster
clusterindex = ismember({EEG.chanlocs.labels}, settings.alphapeak.postelectrodes);
PostSpectro = mean(EEG.specdata(clusterindex,idxlower:idxupper));

%% 1) find the maximimum point
[~ , maxPeak ] = max(PostSpectro);

window = maxPeak-settings.alphapeak.window:maxPeak+settings.alphapeak.window;
if maxPeak - settings.alphapeak.window/freqRes < 1 || maxPeak + settings.alphapeak.window/freqRes > length(freqRange)    
    % this is if the peak is at the borders of the search window
    amplitudeMax = NaN;
    maxPeak = NaN;
else
    amplitudeMax = mean(PostSpectro(window));
end

%% 2) check the first derivative
derivPolarity = diff(PostSpectro) > 0 ; % find the place where it crosses zero coming...
derivPeak = find(diff(derivPolarity) == -1)+1 ; % from positive to negative == peak!
window = derivPeak-settings.alphapeak.window:derivPeak+settings.alphapeak.window;
if isempty(derivPeak)
    derivPeak = NaN;
    amplitudeDer = NaN;
elseif numel(derivPeak) > 1
    derivPeak = NaN;
    amplitudeDer = NaN;
elseif derivPeak - settings.alphapeak.window/freqRes < 1 || derivPeak + settings.alphapeak.window/freqRes > length(freqRange)    
    amplitudeDer = NaN;
    derivPeak = NaN;
else
    amplitudeDer = mean(PostSpectro(window));
end



%% 3) 'gravity method' weighted mean
% uses all electrodes
AllSpectro = mean(EEG.specdata(:,idxlower:idxupper));
%(sum(a(f) x  f))/(Sum a(f)). (Klimesch, 1999)

    gravalphapeak = sum((AllSpectro .* freqRange))./ sum(AllSpectro);

    gravamplitude= interp1(freqRange,AllSpectro,gravalphapeak);

% check if the peak of both methods is equal or not
if maxPeak == derivPeak
    fprintf('\n:::one identifiable peak...\n')
    uniquelyIdentifiable = 1;
elseif isnan(maxPeak) || isnan(derivPeak)
    fprintf('\n:::no identifiable peak...\n')
    uniquelyIdentifiable = 0;
else
    fprintf('\n:::more or no identifiable peak...\n')
    uniquelyIdentifiable = 0;
end

if isnan(maxPeak) 
    out.alphapeakfreqMax = NaN;
    out.alphapeakamplitudeMax = NaN;
    
else
    out.alphapeakfreqMax = freqRange(maxPeak);
    out.alphapeakamplitudeMax = amplitudeMax;
    
end

if isnan(derivPeak)
    out.alphapeakfreqDerivative = NaN;
    out.alphapeakamplitudeDerivative = NaN;
else
    out.alphapeakfreqDerivative = freqRange(derivPeak);
    out.alphapeakamplitudeDerivative = amplitudeDer;
end

out.alphapeakfreqGravity = gravalphapeak;
out.alphapeakamplitudeGravity = gravamplitude;
out.uniquelyIdentifiable = uniquelyIdentifiable;
end


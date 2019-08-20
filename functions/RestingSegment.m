%% this function epochs the data and concatenates it to a continuous EEG file (2D)

function [EEG]  = RestingSegment(EEG,settings)
fprintf('\n:::RestingSegment...\n')

EEG = pop_epoch(EEG,{settings.events},settings.timelimits);

% put it back to continuous EEG format and cleanup the EEG structure
EEG.data = reshape(EEG.data,[EEG.nbchan,EEG.pnts*EEG.trials]);
EEG.pnts = size(EEG.data,2);
EEG.trials = size(EEG.data,3);
EEG.times = 1: 1000/EEG.srate: size(EEG.data,2)*1000/EEG.srate;

end
% calculates the individual frequency bands,based on alpha peak. as
% suggested by: Doppelmayr M, Klimesch W, Pachinger T, Ripper B. Individual
% differences in brain dynamics: important implications for the calculation
% of event-related band power. Biol Cybern 1998;79:49?57: 0.4*IAF?0.6*IAF,
% 0.6*IAF?0.8*IAF, 0.8*IAF?IAF, IAF?1.2*IAF, and 1.2*IAF?25 Hz for theta,
% lower-1-alpha, lower-2-alpha, upper alpha, and beta, respectively
% (Doppelmayr et al. 1998). In addition all possible ratios between
% frequency bands and electrode clusters are computed --> so it is super
% exploratory!!!
% output: Indfbands = individual alpha band power
% This script has not been tested very thoroughly. 

function EEG = IndividualSpectro(EEG,settings)

IAF = EEG.alphaPeak.alphapeakfreqGravity;

Indfbands(1).lowfreqs = 0.4*IAF;
Indfbands(1).highfreqs = 0.6*IAF;
Indfbands(1).name  = 'theta__';
Indfbands(1).absmean  = [];

Indfbands(2).lowfreqs = 0.6*IAF;
Indfbands(2).highfreqs = 0.8*IAF;
Indfbands(2).name  = 'lower_1_alpha__';
Indfbands(2).absmean  = [];

Indfbands(3).lowfreqs = 0.8*IAF;
Indfbands(3).highfreqs = IAF;
Indfbands(3).name  = 'lower_2_alpha__';
Indfbands(3).absmean  = [];

Indfbands(4).lowfreqs = IAF;
Indfbands(4).highfreqs = 1.2 * IAF;
Indfbands(4).name  = 'upper_alpha__';
Indfbands(4).absmean  = [];

Indfbands(5).lowfreqs = 1.2 * IAF;
Indfbands(5).highfreqs = 25;
Indfbands(5).name  = 'beta__';
Indfbands(5).absmean  = [];

Indfbands = computeFbands(EEG.specdata,EEG.freqs,Indfbands,settings.eleclusters,EEG.chanlocs);


 function Indfbands = computeFbands(specdata,freqs,Indfbands,eleclusters,chanlocs)
                
        for f = 1:length(Indfbands)
            % we round to 0.5 because that is our frequency resolution
            ind = find(freqs == round(2*Indfbands(f).lowfreqs)/2) : find(freqs == round(2*Indfbands(f).highfreqs)/2,1);
            Indfbands(f).Roundedlowfreqs = round(2*Indfbands(f).lowfreqs)/2;
            Indfbands(f).Roundedhighfreqs = round(2*Indfbands(f).highfreqs)/2; 
            Indfbands(f).absmean = mean(specdata(:,ind),2);
            
            range = find(freqs == round(2*Indfbands(1).lowfreqs)/2) : find(freqs == round(2*Indfbands(end).highfreqs)/2);
            Indfbands(f).relmean = mean(specdata(:,ind),2)./ mean(specdata(:,range ),2);
            
            % average over electrode clusters
            for  k = 1:length(eleclusters)
                clusterindex = ismember({chanlocs.labels}, eleclusters(k).chans);
                Indfbands(f).elecluster(k).names = eleclusters(k).names;
                Indfbands(f).elecluster(k).absmean = mean(Indfbands(f).absmean(clusterindex));
                Indfbands(f).elecluster(k).relmean = mean(Indfbands(f).relmean(clusterindex));
            end
        end
              
 end

D = [];
% Ratios:
k = 1;
for i1 = 1:length(Indfbands)
    for i2 = 1:length(Indfbands)
        for j1 = 1:length(Indfbands(1).elecluster)
            for j2 = 1:length(Indfbands(1).elecluster)
                if i1 ~= i2 && j1 ~= j2
                D(k).ratio = Indfbands(i1).elecluster(j1).absmean ./ Indfbands(i2).elecluster(j2).absmean;
                D(k).name = [Indfbands(i1).name Indfbands(i2).name Indfbands(1).elecluster(j1).names Indfbands(1).elecluster(j2).names];
                k = k + 1;
                end
            end
        end
    end
end

EEG.Indfbands = Indfbands;
EEG.Ratios = D;

end


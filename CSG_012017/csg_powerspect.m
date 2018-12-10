function varargout = csg_powerspect(varargin)

if nargin < 1
    D = spm_eeg_load;
    eegchan = meegchannels(D);
    epoch = 4;
    flag = 1;
    artefact = [];
elseif nargin == 1    
    if isfield(varargin{1},'Dmeg')
        D = varargin{1}.Dmeg{1};
    end        
    if isfield(varargin{1},'powchan')
        eegchan = varargin{1}.powchan;
    else
        eegchan = meegchannels(D);
    end
    if isfield(varargin{1},'epoch')
        epoch = varargin{1}.epoch;
    else
        epoch = 4;
    end
    if isfield(varargin{1},'plot')
        flag = varargin{1}.plot;
    else 
        flag = 0;
    end
    if isfield(varargin{1},'artefact')
        artefact = varargin{1}.artefact;
    else 
        artefact = [];
    end
end

fprintf(1,'SPECTROGRAM IS BEING COMPUTED \n');
fprintf(1,'==============================\n');


% load parameters
fs = fsample(D);
pow2 = nextpow2(epoch*fs);
nfft = (2^(pow2-1));
nspl = nsamples(D);

% if bad channels - remove them
try 
    badchannels = D.CSG.artefact.badchannels.smallepochs;
    epoch = D.CSG.artefact.badchannels.info.epoch;
catch 
    badchannels = {};
end

% average signal
avg_filt = mean(D(eegchan,:));

% remove bad channels detected from the averaged signal
if ~isempty(badchannels)
    NofE = ceil(nspl/fs/epoch);
    for iep = 1 : NofE
       chan = badchannels{iep};
       if ~isempty(chan)
           goodchan = setdiff(eegchan,chan);
           if isempty(goodchan)
               avg_filt((iep-1)*epoch*fs+1:iep*epoch*fs) = NaN;
           else 
               avg_filt((iep-1)*epoch*fs+1:min(nspl,iep*epoch*fs)) = mean(D(goodchan,(iep-1)*epoch*fs+1:min(nspl,iep*epoch*fs)));
           end
       end 
    end
end

% power by 2 sec epochs (by default)
[S F T P] = spectrogram(avg_filt(:),epoch*fs,[],0.1:fs/nfft:30,fs);
args = {F,T,10*log10(abs(P'))};
if flag
    figure;
    xlbl = 'Frequency (Hz)';
    ylbl = 'Time (sec)';
    surf(args{:},'EdgeColor','none');

    axis xy; axis tight;
    colormap(jet);
    view(90,270);

    ylabel(ylbl);
    xlabel(xlbl);
end

% save power spectrum in the data structure 
D.CSG.spectrogram.info.channels = eegchan;
D.CSG.spectrogram.info.epoch = epoch;
D.CSG.spectrogram.info.artefact = ~isempty(artefact);

D.CSG.spectrogram.tempo = args{2};
D.CSG.spectrogram.frequency = args{1};
D.CSG.spectrogram.power = args{3};

save(D);
varargout{1} = D;


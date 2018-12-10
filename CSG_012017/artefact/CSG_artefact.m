function D = CSG_artefact(cfg)

% detection of artifacts (incoherence) over small epochs after a bad
% channels detection processed over larger time-windows.
% configuration is provided with the structure cfg. :
% °.dataset     : meeg structure 
% °.epoch;      : small time-window length
% °.winsize;    : large time-window length
% °.badchannels; : badchannels detected over large time-windows (after removing channels repaired with the
% interpolation)
% °.channels;   : channels to be analyzed

if ~isempty(cfg)
    D           =   cfg.dataset;
    epoch       =   cfg.epoch;
    winsize     =   cfg.winsize;
    badchannels =   cfg.badchannels;
    channels    =   cfg.channels;
else 
    error('Configuration missing !!!!!!!!!!! \n')
end

% Fix thresholds chosen
tr_spc  =   0.2;     % threshold for spatial coherence: signals have to be coherent within a distance of tr_spc

% load montage
coord    =   coor2D(D,channels);
if size(coord,2)<numel(channels)
    error('coordinates for channels chosen missing !!!!!!!!! \n')
end
Dchan	=   DC_distance(coord');

% data parameters
fs      =   fsample(D);
nspl    =   nsamples(D);
Nchan   =   numel(channels);
Nepo    =   ceil(nspl/fs/epoch);

% initialization
zMatrx      =   NaN(Nepo,Nchan);
zMatrx2     =   NaN(Nepo,Nchan);
A           =   NaN(Nepo,Nchan);
badchan     =   cell(1,Nepo);
Tempo_std   =   std(D(:,1:epoch*fs),[],2);


fprintf(1,'BAD CHANNELS detection over small %d sec-epochs \n', epoch);
fprintf(1,'================================================\n');

for iw  =   2  :    Nepo
    fprintf('.')
    twin    =   (iw-1)*fs*epoch+1:min(nspl,iw*epoch*fs);
    iwin    =   ceil(iw/(winsize/epoch));
    if numel(badchannels{iwin}(:))  >   0.5*Nchan % when more than 50% of channels considered as bad, all channels are 'removed'
        badchan{iw} = channels;
    else 
        for ichan   =   1 : Nchan
        %%  incoherence: kinds of z score computed from small areas for each channel
            chan_around     =   setdiff(channels(Dchan(ichan,:)<tr_spc),[channels(ichan(:)); badchannels{iwin}(:)]);
            if numel(chan_around)>1
                Dmean = mean(D(chan_around,twin));
            else 
                Dmean = D(chan_around,twin);
            end
            zsc = std(D(channels(ichan),twin)-Dmean)/std(Dmean);
            zMatrx(iw,ichan)= zsc;
         %% kinds of zscore over small fix time windows on a same channel
            zsc2 = std(D(channels(ichan),twin));
            zMatrx2(iw,ichan)= zsc2/Tempo_std(channels(ichan));
            if zsc2/Tempo_std(channels(ichan))<3
               Tempo_std(channels(ichan)) = zsc2;
            end  
         %% popping2: compute the maximal slope over small epoch for each channel
            t05 = [1 : 0.5*fs];
            slope = zeros(1,floor(numel(twin)/fs*2));
            for i05 = 1 : floor(numel(twin)/fs*2)
                [Vmax imax] = max(D(channels(ichan),t05),[],2);
                [Vmin imin] = min(D(channels(ichan),t05),[],2);
                slope(i05) = (Vmax-Vmin)/(abs(imax-imin)/fs); 
                t05 = t05 + 0.5*fs;
            end
            A(iw,ichan) = max(slope);
        end
        chantoremove = channels(or(abs(zMatrx(iw,:)>3),abs(zMatrx2(iw,:)>3)));
        if numel(chantoremove)>0.5*Nchan
            badchan{iw} = channels;
        else 
            badchan{iw} = chantoremove;
        end 
    end
end

abn = cell(1,Nchan);
for ichan = 1 : Nchan
    speedchan = A(:,ichan);
    newabn = find(abs(zscore(speedchan))>3);
    count = 0;
    while ~isempty(newabn) && count <= 3
        speedchan(abn{ichan}) = 0;
        newabn =  find(abs(zscore(speedchan))>3);
        abn{ichan} = [abn{ichan}; newabn];
        count = count + 1;
    end
    inepoch = find(diff([0; sort(abn{ichan})])==2);
    abn{ichan} = sort([abn{ichan} ; inepoch]);
    fprintf('.')
end
fprintf('.\n')
% regroup popping ('abn') with badchan...
for ichan = 1 : Nchan
    badepo = abn{ichan};
    for ibe = 1 : numel(badepo)
        badchan{badepo(ibe)} = union(badchan{badepo(ibe)},channels(ichan));
    end
end
for iw = 1 : Nepo
    if numel(badchan{iw})>0.5*Nchan
        badchan{iw} = channels;
    end
end

%% test 
fprintf(1,'---Bad epochs detection --- DONE---\n')
D.CSG.artefact.badchannels.smallepochs = badchan;
D.CSG.artefact.badchannels.info.epoch  = epoch;
save(D);

%%%%
function xf = filterlowhigh(x,frqcut,fs,forder)

flc = frqcut(1)/(fs/2);
fhc = frqcut(2)/(fs/2);
[B,A] = butter(forder,flc,'high');
xf = filtfilt(B,A,x);
[B,A] = butter(forder,fhc,'low');
xf = filtfilt(B,A,xf);

return
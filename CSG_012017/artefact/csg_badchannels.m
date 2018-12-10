
function varargout = csg_badchannels(varargin)

% FORMAT DC_badchannels(args)
% Bad channels detection for whole sleep recordings. 
% Bad channels are detected by scoring window and are either:
%       * Noisy
%      or
%       * Flat
% Badchannels are saved in: 
% * C.CSG.artefact.badchannels.chan_defaillant = flat channels saved in a matrix.
%                           the number of the row correspond to the
%                           EEG channel number, from 1 to the number of EEG channels 
% * C.CSG.artefact.badchannels.chan_incoherent
% the thresholds are available in CRC_get_default('bc')    
%
% INPUT
%       .file   - data file (.mat files)
%__________________________________________________________________
% Copyright (C) 2014 Cyclotron Research Centre

% Written by D. Coppieters 't Wallant, 2014.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

% function to detect badchannels by window of 20s or 30s according to the
% time used for scoring.
% Badchannels are recorded in: 
% * C.CSG.artefact.badchannels.chan_defaillant
% * C.CSG.artefact.badchannels.chan_incoherent
% the thresholds are available in CRC_get_default('bc')
% ----------------------
% *************************************************************************
%                         Loading and parameters values 
% *************************************************************************
if nargin == 1
    D = varargin{1}.dataset;
    winsize =  varargin{1}.winsize;
    channels = varargin{1}.channels;
else 
    D = spm_eeg_load;
    winsize = 30;
    channels    =   meegchannels(D);
end


fs = fsample(D);   
nspl    =   nsamples(D);
Time    =   nspl/fs;
NofW    =   ceil(Time / winsize);
nbr_sc_epoch    =   winsize*fs;

% *************************************************************************
%                                 Threshold 
% *************************************************************************
clear global crc_def;
tr = crc_get_defaults('qc.bc');
% Obvious channels --------------------------------------------
tr_n  = tr.n;  % noisy channel threshold
tr_f1 = tr.f1; % 1st threshold of flat channel
% Noisy channels ----------------------------------------------
tr_r  = tr.r; % ratio of deviation
% Flat channels -----------------------------------------------
tr_f2 = tr.f2; % 2nd threshold of flat channel
tr_tf = tr.tf; % duration of flat channel
tr_ampl = tr.ampl; % amplitude depending on standard deviation

% initialization
chan_def = cell(NofW,1);
chan_incoh = cell(NofW,1);

% *************************************************************************
%                    Obviously bad channels detection
% *************************************************************************
h = waitbar(0,'badchannels detection');  

fprintf(1,'BAD CHANNELS detection over large %d sec-epochs \n', winsize);
fprintf(1,'================================================\n');
for w = 1 : NofW
    window = D(channels,(w-1)*nbr_sc_epoch+1 : min(w*nbr_sc_epoch,nspl));
       % --- obvious noisy channel
       eeg_an    =   std(window,[],2); 
       if any(eeg_an >= tr_n)
           idc = find(eeg_an >= tr_n);
           for i = 1 : length(idc)
                chan_incoh{w} = channels(idc);
           end

       % --- obvious flat channel
       elseif any(eeg_an <= tr_f1) 
           idc = find(eeg_an <= tr_f1);
           for i = 1 : length(idc)
               chan_def{w} = channels(idc);
           end
       end
    String  =  ['Progress of bad channels detection 1/2 : ' num2str(w/NofW*100) ' %'];
    waitbar((w/NofW),h,String);
end
% *************************************************************************
%                Finer detection for noisy and flat channels
% *************************************************************************
for w = 1 : NofW
    cki = chan_incoh{w};
    ckn = chan_def{w};
    good_chan = setdiff(channels,[cki;ckn]); % index of obvious bad channels = [cki and ckn] so we keep index of potential good channels
    if ~isempty(good_chan)
       window     =   D(good_chan,(w-1)*nbr_sc_epoch+1 : min(w*nbr_sc_epoch,nspl));
       st_5epo    =   std(window,[],2); % standard deviation along each channel
       for c = 1:numel(good_chan)
            w_c = window(c,:);         
% ***** flat channel ***** (Devuyst)
            if std(window(c,:))<tr_f2
                def = abs(w_c);
                v = find(def <= tr_ampl*st_5epo(c)); 
                diff_def = diff(v);
                fin = unique([find(diff_def~=1) length(diff_def)]);
                deb = [1 fin(1:end-1)+1];
                duration_def = (fin - deb)/fs;
                if any(duration_def >= tr_tf)
                    chan_def{w} = union(chan_def{w},good_chan(c));  %artefact de défaillance des électrodes
                end  
            end
% ***** noisy channels *****
            other     =     setdiff(1:numel(good_chan),c); %all good channels except the one we are analyzing
            eeg_oth   =  	mean(window(other,:));
            eeg_an    =     window(c,:); 
            intera_std =    abs(std(eeg_an)/std(eeg_oth)); 
            if (intera_std >= tr_r)
               chan_incoh{w} = union(chan_incoh{w},good_chan(c));  %artefact de défaillance des électrodes
            end
       end
    end
    String  =  ['Progress of bad channels detection 2/2: ' num2str(w/NofW*100) ' %'];
    waitbar((w/NofW),h,String);
end
close(h);
D.CSG.artefact.badchannels.chan_defaillant = chan_def;
D.CSG.artefact.badchannels.chan_incoherent = chan_incoh;
D.CSG.artefact.badchannels.info.winsize = winsize;
D.CSG.artefact.badchannels.info.chan = chanlabels(D,channels);
save(D);
varargout{1} = D;
fprintf(1,'---Bad channels detection DONE---\n')
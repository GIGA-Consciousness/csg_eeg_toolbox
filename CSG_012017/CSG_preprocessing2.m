function Dprep = CSG_preprocessing2D

% Complete preprocessing
clear;
clc;

%% Convert to the SPM format 
% % S.dataset = spm_select(1,'any','Select raw data file');
% % S.outfile = S.dataset;
% % D = spm_eeg_convert(S);

D = spm_eeg_load;
% parameters 
fs      =   fsample(D);
nspl    =   nsamples(D);

%% Clone data ---------
% newfilename  =  fullfile(path(D),['P_' fname(D)]);
% Dnew    =   clone(D,newfilename);

%% Filter 
% -------
% 
eegchan = meegchannels(D);
% h = waitbar(0,'EEG channels filtering'); 
% for ichan = 1 : numel(eegchan)
%     waitbar((ichan/numel(eegchan)),h,sprintf('%d/%d filtered EEG channels', ichan, numel(eegchan)));
%     Dnew(eegchan(ichan),:) = filterlowhigh(D(eegchan(ichan),:),[0.1 30],fs,3);
% end
% eogchan = eogchannels(D);
% if ~isempty(eogchan) 
%     for ichan = 1 : numel(eogchan)
%         waitbar((ichan/numel(eogchan)),h,sprintf('%d/%d filtered EOG channels', ichan, numel(eogchan)));
%         Dnew(eogchan(ichan),:) = filterlowhigh(D(eogchan(ichan),:),[0.1 5],fs,3);
%     end
% end
% emgchan = emgchannels(D);
% if ~isempty(emgchan) 
%     for ichan = 1 : numel(emgchan)
%         waitbar((ichan/numel(emgchan)),h,sprintf('%d/%d filtered EMG channels', ichan, numel(emgchan)));
%         Dnew(emgchan(ichan),:) = filterlowhigh(D(emgchan(ichan),:),[10 125],fs,3);
%     end
% end
% otherchan = setdiff(1:size(D,1),[eegchan eogchan emgchan]);
% if ~isempty(otherchan) 
%     for ichan = 1 : numel(otherchan)
%         waitbar((ichan/numel(otherchan)),h,sprintf('%d/%d filtered OTHER channels', ichan, numel(otherchan)));
%         Dnew(otherchan(ichan),:) = filterlowhigh(D(otherchan(ichan),:),[0.1 30],fs,3);
%     end
% end
% close(h);
% save(Dnew);

%% BAD CHANNELS over large fix time-windows
% -----------------------------------------

cfg = [];
cfg.dataset = D;
cfg.winsize = 20;
cfg.channels = eegchan;
% 
% Dbadc = csg_badchannels(cfg);
% 
% %% Interpolation of bad channels detected
% % ---------------------------------------
NofW    =   ceil(nspl/(fs*cfg.winsize));
% for ibc = 1 : NofW
%     cfg.badchannels{ibc} = union(Dbadc.CSG.artefact.badchannels.chan_defaillant{ibc},Dbadc.CSG.artefact.badchannels.chan_incoherent{ibc});
% end
% cfg.dataset      =   Dbadc;
% [Drbadc cleaned] =   csg_interpol(cfg);

% keep in memory channels repaired and remove them from badchannels
% Drbadc.CSG.artefact.badchannels.largecleaned = cleaned;
% for ibc = 1 : NofW
%     Drbadc.CSG.artefact.badchannels.chan_defaillant{ibc} = setdiff(Drbadc.CSG.artefact.badchannels.chan_defaillant{ibc},Drbadc.CSG.artefact.badchannels.largecleaned{ibc});
%     Drbadc.CSG.artefact.badchannels.chan_incoherent{ibc} = setdiff(Drbadc.CSG.artefact.badchannels.chan_incoherent{ibc},Drbadc.CSG.artefact.badchannels.largecleaned{ibc});
% end
% save(Drbadc);
% cfg.dataset  =  Drbadc;

%% BAD CHANNELS detection over small epochs
% -----------------------------------------
cfg.epoch   =   4;
% % epowin      =   cfg.winsize/cfg.epoch;
for ibc = 1 : NofW
    cfg.badchannels{ibc} = union(D.CSG.artefact.badchannels.chan_defaillant{ibc},D.CSG.artefact.badchannels.chan_incoherent{ibc});
end
Dprep           =   CSG_artefact(cfg);
% % cfg.dataset     =   Dprep;
% % cfg.badchannels =   cfg.dataset.CSG.artefact.badchannels.smallepochs;
% % cfg.winsize     =   cfg.dataset.CSG.artefact.badchannels.info.epoch;

% % %% interpolation of bad channels detected over small epochs
% % % [Dprep cleaned] = csg_interpol(cfg);
% % 
% % % keep in memory channels repaired and remove them from badchannels
% % % NofE    =   ceil(nspl/(fs*cfg.winsize));
% % % Dprep.CSG.artefact.badchannels.smallcleaned = cleaned;
% % % for ibc = 1 : NofE
% % %     Dprep.CSG.artefact.badchannels.smallepochs{ibc} = setdiff(Dprep.CSG.artefact.badchannels.smallepochs{ibc},Dprep.CSG.artefact.badchannels.smallcleaned{ibc});
% % % end

% add bad channels found over large time windows to those detected in small
% time windows
% % badchannels = cell(1,NofW);
% % for ibc = 1 : NofW    
% %     badchannels{ibc} = union(Dprep.CSG.artefact.badchannels.chan_defaillant{ibc},Dprep.CSG.artefact.badchannels.chan_incoherent{ibc});
% %     epochs = (ibc-1)*epowin + 1 : min(NofE,ibc*epowin);
% %     for ie = epochs
% %         Dprep.CSG.artefact.badchannels.smallepochs{ie} = union(Dprep.CSG.artefact.badchannels.smallepochs{ie},badchannels{ibc});
% %     end
% % end

save(Dprep);


%%%%%%%%%%%%%%%%%%%%
function xf = filterlowhigh(x,frqcut,fs,forder)

flc = frqcut(1)/(fs/2);
fhc = frqcut(2)/(fs/2);
[B,A] = butter(forder,flc,'high');
xf = filtfilt(B,A,x);
[B,A] = butter(forder,fhc,'low');
xf = filtfilt(B,A,xf);

return


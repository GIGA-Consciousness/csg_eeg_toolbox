% Complete preprocessing
clear;
file = spm_select(1,'any','Select raw data file');

% 1) ft-preprocessing to creater the header and preprocess data

%%% filtering configuration %%%
% -----------------------------
D = crc_eeg_load(file);
newfilename = fullfile(path(D),['P' fname(D)]);
Dnew = clone(D,newfilename);
data = spm2fieldtrip(D);

cfg.lpfilter = 'yes';
cfg.hpfilter = 'yes';
cfg.lpfreq = 30;    % Low pass filter
cfg.hpfreq = 0.5;   % High pass filter
%%%
[pth,dataname,ext] = fileparts(fname(D));
cfg.datafile = fnamedat(D);
cfg.headerfile = fullfile(path(D),['P' fname(D)]);
cfg.outputfile = fullfile(path(D),['P' fname(D)]);
cfg.channel = chanlabels(D);
newfile = ft_preprocessing(cfg,data);
Dnew(:,:,1) =  newfile.trial{1};
save(Dnew);

% 
% cfg2.method = 'spline';
% cfg2.elecfile = 'C:\Users\Doro\Github\EEGtoolswc\branches\brancheDC\FASST\CSG\Test\elec_EGI256.mat';
% [test] = ft_scalpcurrentdensity(cfg2, newfile)


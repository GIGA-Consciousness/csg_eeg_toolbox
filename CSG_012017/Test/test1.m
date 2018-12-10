function D = csg_interpol(varargin)

% csg_interpol interpolate bad channels detected from the csg_badchannels detection method (for each 30s window).
% if varargin is not empty, it is a struct containing:
% ° .filename*   : filename of a meeg file
% ° .method      : 'spline' or 'average' interpolation method
% ° .lambda      : smoothing (?)
% ° .order       : m-constant for the spline interpolation
%  This code needs functions fromt the fieldtrip toolbox (i.e. ft_channelrepair) available in SPM12

if nargin>0
    D = spm_eeg_load(varargin{1}.filename);
else 
    D = spm_eeg_load;
end

% to load the sensors positions of the geodesic sensors with 256 channels
load('C:\Users\Doro\Github\EEGtoolswc\branches\brancheDC\FASST\CSG\Test\elec_EGI256.mat')

% parameters from the data file loaded
fs      =   fsample(D);        % sampling frequency
Lw      =   30;                % fix time window used for bad channels detection
nspl    =   nsamples(D)        % number of samples
NofW    =   ceil(nspl/(fs*Lw));% number of 30s window
fcl     =   0.5;               % low cut frequency
fch     =   30;                % high cut frequency
forder  =   3;                 % filter order

badchannels = DC_union(D.CSG.artefact.badchannels.chan_defaillant,D.CSG.artefact.badchannels.chan_incoherent);
data = spm2fieldtrip(D); % transform meeg object from spm into a raw data of fieldtrip
newdata = data; % new raw data with interpolated values 

% data parameters
fs      = fsample(D);
nspl    = nsamples(D);
Nchan   = size(D,1);
Nepo = ceil(nspl/fs/epoch);

% filtering
fD = zeros(size(D));
fprintf(1,'1) Filtering data \n');
for ichan = 1 : Nchan
    fD(ichan,:) = filterlowhigh(D(ichan,:),[fcl fch],fs,forder);
end

% check the configuration chosen and put default value if not available
% values 
cfg.elec    =   elec; % channels position and labels given from the elec_EGI256.mat file according to the fieldtrip toolbox
cfg.method  =   ft_getopt(cfg, 'method','spline');
cfg.lambda  =   ft_getopt(cfg, 'lambda',[]); % subfunction will handle this
cfg.order   =   ft_getopt(cfg, 'order',[]); % subfunction will handle this

% loop to search bad channels in each 30s-window and interpolate bad
% channels from the method chosen

for iw = 1 : NofW

    cfg.badchannel = chanlabels(D,badchannels{iw});
    newdata.trial{1} = fD(:,(n-1)*Lw*fs:n*fs*Lw);
data_test.time{1} = data.time{1}(:,(n-1)*Lw*fs:n*fs*Lw);
data_test2 = ft_channelrepair(cfg,data_test);



%%%%
function xf = filterlowhigh(x,frqcut,fs,forder)

flc = frqcut(1)/(fs/2);
fhc = frqcut(2)/(fs/2);
[B,A] = butter(forder,flc,'high');
xf = filtfilt(B,A,x);
[B,A] = butter(forder,fhc,'low');
xf = filtfilt(B,A,xf);

return

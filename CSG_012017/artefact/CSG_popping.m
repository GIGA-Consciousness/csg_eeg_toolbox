function varargout = csg_popping(varargin)

% FORMAT csg_popping(args)
% It detects rapid transitions also called poppin artifacts over all good EEG channels  
% The rapid transitions are evaluated by 0.5-s epochs.
%        popping artifacts are saved in C.CSG.artefact.popping
% INPUT
%       .file   - data file (.mat files)
% 
% NB: it has to be processed after DC_preprocessing, DC_badchannels and DC_pwr
% detection) in every 1-s epoch
% _________________________________________________________________
% Copyright (C) 2014 Cyclotron Research Centre

% Written by D. Coppieters 't Wallant, 2016.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$
% ----------------------

if nargin == 1
    D = spm_eeg_load(varargin{1});
else 
    D = spm_eeg_load;
end

winsize =   30;
fs      =   fsample(D);
nspl    =   nsamples(D);
Time    =   ceil(nspl/fs);
channels = meegchannels(D);

% Initialization 
ae = zeros(1,2*Time);
te = zeros(1,2*Time);
Re = zeros(1,2*Time);
for epoch = str
    win = ceil(epoch/winsize);
    [ch c v] = find(badchannels==win);
    gd = setdiff(channels,ch);
    if ~isempty(gd)
        for it = (epoch-1)*2+1 : min(ceil(2*nspl/fs),epoch*2)
            treeg = C(gd,((it-1)*fs/2+1: min(nspl,it*fs/2)));
            [vmax lmax] = max(treeg,[],2);
            [vmin lmin] = min(treeg,[],2);
            [Re(it) c] = max((vmax-vmin)./(abs(lmax-lmin)/fs));
            ae(it) = vmax(c)-vmin(c);
            te(it) = (abs(lmax(c)-lmin(c))/fs);
            String  =  ['Progress for popping artefact on all EEG channels... : ' num2str(100*epoch/(Time)) ' %'];
            waitbar(epoch/(Time),h,String);
        end   
    end
end  
popping = unique(ceil(find(and(Re>3e3,or(ae<120,te<0.03)))/2));
close(h) %enregistrement------------------------------------------------------------
C.CRC.DC.shortartf.popping = popping;
save(C);
fprintf('* popping artifacts detected per 1-s epoch: %d \n',length(popping))

function varargout = csg_reref(varargin)
%
% Easy rereferencing of EGI 256 channels to mean mastoids: 
% Left Mastoid : E91
% Right Mastoid : E216
%
% WARNING: you must launch "spm eeg" before using this function as it
% relies on SPM functionalities (referencing function and fieldtrip tools).
%__________________________________________________________________
% Copyright (C) 2016 Coma Science Group

% Written by D. Coppieters 
% Coma Science Group, University of Liege, Belgium
% $Id: $

% REF2 = mean(E91,E216);
% X - REF2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Some definitions:
%------------------
origdir = 'C:\Users\Doro\Github\EEGtoolswc\branches\brancheDC\FASST\CSG';
montage_file = 'Montage_CSG_256_avgall.mat';

% montage options
%----------------
S.keepothers        = 1; 
S.updatehistory     = 1;
S.onlineopt         = 0; % apply montage and write out data

% getting things started
%-----------------------
cd(origdir)
load(montage_file)
if nargin == 1
    files = varargin{1}.file;
else 
    files = spm_select(Inf,'mat','Select subject''s files');
end
% Looping through all the subjects
%---------------------------------
for iP = 1:size(files,1)
    D = spm_eeg_load(deblank(files(iP,:)));
    S.D                 = D;
    S.montage.tra       = montage.tra;
    S.montage.labelorg  = montage.labelorg;
    S.montage.labelnew  = montage.labelnew;
    if isfield(montage,'name');
        S.montage.name  = montage.name;
    end

    mD = spm_eeg_montage(S);
    D = chantype(mD,[91 216],'REF');
    save(D);
end
cd(origdir)
varargout{1} = D;
fprintf(1,'---Rereferencing DONE---\n')

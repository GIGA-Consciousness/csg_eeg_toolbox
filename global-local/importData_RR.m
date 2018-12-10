%% Clear the workspace
%ft_defaults                            % only at the beginning of each Matlab session
clear all
clc
SUBJECTpath = spm_select(1,'.mat','Select the EEG file'); % Subject folder
D = spm_eeg_load(SUBJECTpath);
cd(path(D));
    
%% ----------------------------------------------------------------------------
% Define trials
blocks = {1:8};  % Which blocks are in each dataset? here the 8 blocks were recorded in a unique dataset %%{1:4, 5:8};                

% Load the stimulation sequence and recreate the conditions   %load([PATH, SUBJECT, 'GlobLoc_stimulation_SEQUENCE'], 'SEQUENCE')
GL_stimulation = ls('GlobLoc_stimulation_SEQUENCE*');
if size(GL_stimulation,1)>1
    fprintf('Only one ''GlobLoc_stimulation_SEQUENCE'' file should be present in this folder')
    GL_stimulation = spm_select(1,'.mat','Select the ''GlobLoc_stimulation_SEQUENCE'' file', [],pwd,'GlobLoc_stimulation_SEQUENCE*');
end
load(GL_stimulation,'SEQUENCE')

CHORDS = {};
for ds = 1:length(blocks)
    
    CHORDS{ds} = [];
    
    for i = blocks{ds}
    
        standard = SEQUENCE{1,i}(1,:);
    
        if isequal(standard, [1 1 1 1 1])
           block = 1;
       elseif isequal(standard, [2 2 2 2 2])
           block = 2;
       elseif isequal(standard, [1 1 1 1 2])
           block = 3;
       elseif isequal(standard, [2 2 2 2 1])
           block = 4;     
        end
   
        thc = repmat(block,size(SEQUENCE{1,i},1),5);
        thc(SEQUENCE{1,i}(:,5) ~= standard(5),5) = block*-1;
        thc = reshape(thc.',1,[]);
    
        % for the 5th triggers afterwards
        trg = repmat([0 0 0 0 1],size(SEQUENCE{1,i},1),1);
        trg = reshape(trg.',1,[]);
    
        CHORDS{ds} = [CHORDS{ds}; [thc' trg']];
    end
end

%% Identify triggers and check them
cfg = [];
cfg.trialdef.eventtype  = 'falling';     %'Response';
cfg.trialdef.eventvalue = 2;         %R128';
cfg.continuous          = 'yes';   
% cfg.trialdef.prestim    = 0.9; % in seconds
% cfg.trialdef.poststim   = 0.8; % in seconds

cfgs = {}; diffs = {};
for ds = 1:length(blocks)
    cfg.dataset = SUBJECTpath;%[PATH, SUBJECT, DATSETS{ds}, '.mat']; %'.vhdr'];
    cfgs{ds}    = ft_definetrial(cfg);
    D = spm_eeg_load(cfg.dataset);
    fs = fsample(D);
    diffs{ds} = diff(cfgs{ds}.trl(:,1)) /fs;   
end
%% Look at the diffs and decide whether some triggers at the beginning do not belong to the paradigm

% should be around a sequence of aprox. 0.15 0.15 0.15 0.15 >1.0

% 1st is the test trigger?
% cfgs{...}.trl(1,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate the missing triggers based on psychtoolbox output

% This is the psychtoolbox output for each stimulus
GL_critical_event = ls('GlobLoc_critical_events*');%load([PATH, SUBJECT, 'logs/GlobLoc_critical_events'])
if size(GL_critical_event,1)>1
    GL_critical_event = spm_select(1,'.mat','Select the ''GlobLoc_critical_events'' file', [], pwd,'GlobLoc_critical_events*');
end
load(GL_critical_event)
count = 1;
missTrig = {};
events   = cfgs{ds}.trl(:,1);
for ds = 1:length(blocks)
    
    tzz = tz(count:count+size(CHORDS{ds},1)-1); % diff(tz);
    
    TZ_line = zeros(1,round((tz(end)-tz(1))*1450)+100); 
    TZ_line(round((tz-tz(1))*1449.3) +1 ) = 2;            % -.7 correction of sampling rate : fs adjusted
    TZ = find(TZ_line)' + cfgs{1}.trl(1,1)-1;

    events   = cfgs{ds}.trl(:,1);
    
    missTrig{ds} = [];
    for i = 1:size(TZ,1)
        if i > numel(events)
            missTrig{ds} = [missTrig{ds}; i];
            vai2ins = TZ(i,1) + events(i-1,1) - TZ(i-1,1); % estimate the sample number
            events = [events(1:i-1);   vai2ins];
            cfgs{ds}.trl = [cfgs{ds}.trl(1:i-1,:);   [vai2ins vai2fins+(cfgs{ds}.trl(1,2)-cfgs{ds}.trl(1,1)) ...
                cfgs{ds}.trl(1,3) cfgs{ds}.trl(1,4)]];
        end
        if TZ(i,1)>events(i,1)-0.15*fs/2 && TZ(i,1)<events(i,1)+0.15*fs/2 % if the sample number is about the same % 50 for 2500 Hz
            continue % Do nothing
        else
            missTrig{ds} = [missTrig{ds}; i];
            vai2ins = TZ(i,1) + events(i-1,1) - TZ(i-1,1); % estimate the sample number
            events = [events(1:i-1);   vai2ins;    events(i:end)];
            cfgs{ds}.trl = [cfgs{ds}.trl(1:i-1,:);   [vai2ins vai2ins+(cfgs{ds}.trl(1,2)-cfgs{ds}.trl(1,1)) ...
              cfgs{ds}.trl(1,3) cfgs{ds}.trl(1,4)];    cfgs{ds}.trl(i:end,:)];
        end
    end
    count = count + size(CHORDS{ds},1);
end

%% Define only the 5th tone triggers as trials
for ds = 1:length(blocks)
    cfgs{ds}.trl = cfgs{ds}.trl(CHORDS{ds}(:,2)==1,:); 
end

[pth, name] = fileparts(GL_stimulation);
namefile = name(numel('GlobLoc_stimulation_SEQUENCE')+1:end);
save(['GlobLoc_trl_', namefile],'CHORDS')

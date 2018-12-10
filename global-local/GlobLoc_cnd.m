%ft_defaults
clear all
clc
%%
PATH    = 'C:\Users\Doro\Documents\Projects\GlobalLocal\';   % Project folder
SUBJECT = 'Patient_0216NC\';                                 % Subject folder

load([PATH, SUBJECT, 'GlobLoc'])
load([PATH, SUBJECT, 'GlobLoc_trl'], 'fake_trl_raw', 'finals', 'CHORDS')

CHORDS = cell2mat({cat(1, CHORDS{:})});
CHORDS = CHORDS(CHORDS(:,2)==1,1); 

[c,ia,ib] = intersect(fake_trl_raw(:,1),finals(:,1)-100);       % ia == millised read leiduvad mõlemas vektoris
CNDs = CHORDS(ia);

%% Conditions and averaging 

cfg = [];
cfg.detrend        = 'no';                        % remove linear trends from the data
cfg.demean         = 'yes';                       % baseline correction
cfg.baselinewindow = [-0.1 0.0];

cfg.trials     = find(CNDs == 1 | CNDs == 2);
LSGS      = ft_preprocessing(cfg, good_data);	% does the preprocessing with the above options
cfg.trials     = find(CNDs == 3 | CNDs == 4);
LDGS      = ft_preprocessing(cfg, good_data);
cfg.trials     = find(CNDs == -1 | CNDs == -2);
LDGD        = ft_preprocessing(cfg, good_data);	% does the preprocessing with the above options
cfg.trials     = find(CNDs == -3 | CNDs == -4);
LSGD        = ft_preprocessing(cfg, good_data);

cfg.trials     = find(CNDs == 1 | CNDs == 2 | CNDs == -3 | CNDs == -4);
LS = ft_preprocessing(cfg, good_data);	% does the preprocessing with the above options
cfg.trials     = find(CNDs == 3 | CNDs == 4 | CNDs == -1 | CNDs == -2);
LD = ft_preprocessing(cfg, good_data);
cfg.trials     = find(CNDs > 0);
GS        = ft_preprocessing(cfg, good_data);
cfg.trials     = find(CNDs < 0);
GD        = ft_preprocessing(cfg, good_data);	% does the preprocessing with the above options
    

cfg = [];
cfg.keeptrials = 'no';
LSGS_avg = ft_timelockanalysis(cfg, LSGS);
LDGS_avg = ft_timelockanalysis(cfg, LDGS);
LSGD_avg = ft_timelockanalysis(cfg, LSGD);
LDGD_avg = ft_timelockanalysis(cfg, LDGD);

cfg = [];
cfg.keeptrials = 'no';
LS_avg = ft_timelockanalysis(cfg, LS);
LD_avg = ft_timelockanalysis(cfg, LD);
GS_avg = ft_timelockanalysis(cfg, GS);
GD_avg = ft_timelockanalysis(cfg, GD);

clear('good_data')

%% Difference waves

LdiffGS = LSGS_avg;
LdiffGS.avg = LDGS_avg.avg - LSGS_avg.avg;
LdiffGD = LSGD_avg;
LdiffGD.avg = LDGD_avg.avg - LSGD_avg.avg;

LSGdiff = LSGS_avg;
LSGdiff.avg = LSGD_avg.avg - LSGS_avg.avg;
LDGdiff = LDGS_avg;
LDGdiff.avg = LDGD_avg.avg - LDGS_avg.avg;


Ldiff = LS_avg;
Ldiff.avg = LD_avg.avg - LS_avg.avg;
Gdiff = GS_avg;
Gdiff.avg = GD_avg.avg - GS_avg.avg;
%% Initial plotting 

cfg = [];
cfg.layout = 'M1_XYZ3_new.sfp';
figure
ft_multiplotER(cfg, LSGS_avg, LDGS_avg, LSGD_avg, LDGD_avg);
figure
ft_multiplotER(cfg, LdiffGS, LdiffGD, LSGdiff, LDGdiff);

figure
ft_multiplotER(cfg, LS_avg, LD_avg);
figure
ft_multiplotER(cfg, GS_avg, GD_avg);
figure
ft_multiplotER(cfg, Ldiff, Gdiff);

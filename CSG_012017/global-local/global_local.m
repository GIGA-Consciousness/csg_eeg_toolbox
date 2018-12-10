%% The global-local paradigm

% Created by Renate Rutiku; last tested December 9th, 2016
% The GenerateTone function is from another source, but cannot find it anymore

% Most parameters are replicated from Faugeras et al. (2012)

% Stimuli: 
% - chords of 3 sinusoids (either 350, 700, and 1400 Hz; or 500, 1000, and 2000 Hz)
% - second and third partials are of 1/2 and 1/4 intensity, respectively
% - duration 50 ms (including 7 ms rise and fall times)
% - presented in series of 5 with 150 ms SOA between sounds
% - series are separated by a variable interval of 1350 to 1650 ms (50-ms steps)

% Conditions:
% - the standard is either a series of the same 5 sounds (AAAAA or BBBBB) or of the final sound swapped (either AAAAB or BBBBA)
% - the deviants' 5th sound is opposite to their respective standard's 5th sound

% Blocks:
% - each block consists of 80% standards and 20% deviants
% - the number of deviants varies randomly between 22 and 30
% - all deviants are separated by at least 2 standards
% - each block starts with 20 –30 standards

% Sequence:
% - each of the 4 possible block types is presented twice
% - the blocks are presented in in a fixed order (two runs of AAAAA, BBBBB, AAAAB, BBBBA global standards)

%%

clear; clc;                              % clear MATLAB environment

mkdir('logs')                            % create a log folder for the files

timestamp = datestr(now);                % date and time to append to the file names
timestamp = regexprep(timestamp, '-','');
timestamp = regexprep(timestamp, ' ','');
timestamp = regexprep(timestamp, ':','');

diary(['logs/GlobLoc_diary_', timestamp])

%% generate the sounds

sf         = 44100;                      % sampling frequency (Hz); make sure that's what the audio card is running
sDuration  = 50;                         % sound duration (ms)
Freqs_low  = [350,700,1400];             % low tone frequencies (Hz)
Freqs_high = [500,1000,2000];            % high tone frequencies (Hz)
sAmps      = [1, 1/2, 1/4];              % normalized chord amplitude modulations (dB)

% Cosine ramp; from http://www.h6.dion.ne.jp/~fff/old/technique/auditory/matlab.html 
dr = 0.007;                              % 7 ms ramps
nr = floor(sf * dr);
CSramp = sin(linspace(0, pi/2, nr));
CSramp = [CSramp, ones(1, floor(sf * sDuration/1000) - nr * 2), fliplr(CSramp)];


% The low tone
lTone = GenerateTone(sf, sDuration, Freqs_low, sAmps); %
lTone = lTone' .* CSramp;                % ramp sound       %plot(lTone)

% The high tone
hTone = GenerateTone(sf, sDuration, Freqs_high, sAmps);
hTone = hTone' .* CSramp;                %plot(hTone)

stims = {};
stims{1,1} = lTone;                      % 1. index == low tone (A)
stims{1,2} = hTone;                      % 2. index == high tone (B)

%% Create the stimulation sequence

Blocks = sort(repmat(1:4,1,2));

SEQUENCE = {};

for bl = 1:numel(Blocks)
    
   SEQUENCE{bl} = [];
   
   % the standards and deviants for each block type; 
   % numbers will reference the "stims" object later on
   if Blocks(bl) == 1
       standard = [1 1 1 1 1];
       deviant = [1 1 1 1 2];
   elseif Blocks(bl) == 2
       standard = [2 2 2 2 2];
       deviant = [2 2 2 2 1];
   elseif Blocks(bl) == 3
       standard = [1 1 1 1 2];
       deviant = [1 1 1 1 1];
   elseif Blocks(bl) == 4
       standard = [2 2 2 2 1];
       deviant = [2 2 2 2 2];       
   end
   
   % introduce the standard before the first deviant
   habit = randsample(20:30,1);
   SEQUENCE{bl} = [SEQUENCE{bl}, repmat(standard, habit,1)];
   
   % how many deviants and 4 times more standards == 20/80
   nrDeviants = randsample(22:30,1);
   nrStandard = 4*nrDeviants;
   
   % order of stimuli
   b2rand = [ones(nrStandard,1); ones(nrDeviants,1)*2];
   b2rand = b2rand(randperm(numel(b2rand)));
   % make sure there are at least 2 standards between consecutive deviants
   while ~isempty(findstr(b2rand',[2 2])) || ~isempty(findstr(b2rand',[2 1 2]))  
       TT= findstr(b2rand',[2 2]);
       if ~isempty(TT)
            b2rand = [b2rand(1:TT(1)); [1;1]; b2rand(TT(1)+1:end)];
            b2rand(randsample(findstr(b2rand',ones(1,7)),2)) = [];
       end
       TOT= findstr(b2rand',[2 1 2]);
       if ~isempty(TOT)
            b2rand = [b2rand(1:TOT(1)); 1; b2rand(TOT(1)+1:end)];
            b2rand(randsample(findstr(b2rand',ones(1,6)),1)) = [];
       end
   end
   
   % add to the SEQUENCE
   for i = 1:numel(b2rand)
       if b2rand(i) == 1
           SEQUENCE{bl} = [SEQUENCE{bl}; standard];
       elseif b2rand(i) == 2
           SEQUENCE{bl} = [SEQUENCE{bl}; deviant];
       end
   end
   
end
% sum(SEQUENCE{1,2}(:,5) ~= SEQUENCE{1,2}(1,1)) % sanity checks

save(['logs/GlobLoc_stimulation_SEQUENCE_', timestamp], 'SEQUENCE')

%% Run the sequence 

% PsychPortAudio('Close')
% Running on PTB-3? Abort otherwise.
AssertOpenGL;
% Perform basic initialization of the sound driver:
InitializePsychSound(1);

% Open the default audio device [], with default mode [] (==Only playback),
% and a required latencyclass of 1 == friendly low-latency mode
% a frequency of sf and 1 sound channel
% This returns a handle to the audio device:
pahandle = PsychPortAudio('Open', [], [], 1, sf, 1);


pause(2.0) % short delay before the sequence begins 

tz = zeros(1,9000);
count = 1;

tic;
for bl = 1:numel(Blocks)
    for trl = 1:size(SEQUENCE{bl},1)

        for i = 1:5
            PsychPortAudio('FillBuffer', pahandle, stims{1,SEQUENCE{bl}(trl,i)});
            tz(count) = PsychPortAudio('Start', pahandle, [], 0, 1); count = count+1;
%             UseTriggerBox('Send',1,0,0,0,0,0,0,0);
%             UseTriggerBox('Send',0,0,0,0,0,0,0,0);
            java.lang.Thread.sleep(140);  % in ms
        end
        
        java.lang.Thread.sleep(randsample(650:50:950,1));

    end
    
    if bl == 1
        [secs, keyCode, deltaSecs] = KbStrokeWait; % After the 4th block wait for a keypress to continue
    else
        pause(5.0)                                 % short delay before the next block begins
    end
end
toc;

PsychPortAudio('Close')
diary off

tz = tz(1:find(tz,1,'last'));                      % clear zeros
% SAVE!!!!
save(['logs/GlobLoc_critical_events', timestamp], 'tz')

%%
save(['logs/GlobLoc_backup_dump_', timestamp])

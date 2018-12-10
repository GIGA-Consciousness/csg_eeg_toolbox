
First things first:
Install psychtoolbox: http://psychtoolbox.org/download/
Change paths in global_local.m according to your system (e.g., windoze uses "\" instead of "/")


Included in the folder:
1) global_local.m ... Runs the paradigm
2) GenerateTone.m ... A function for creating tones
3) randsample.m   ... a MATLAB function; in case it's missing from the installed version of MATLAB


After running global_local.m you should have 4 files in the log folder:
1) a diary of the command line output
2) the stimulation sequence, i.e., a complete specification of the sounds played
3) timestamps for all sound onsets returned by psychtoolbox
4) a backup dump of the entire workspace after the script has finished

######################################################################
Important:
* The sampling frequency (variable sf in the code) should be the same that the audio card runs. 44100 Hz is the most typical value, but please check anyway. 
* If you cancel the sequence before it finishes and still want to save the timing data created so far you have to run the saving functions at the end of the code manually. 
* If you have cancelled a running sequence and want to start from the beginning run PsychPortAudio('Close') first. Otherwise calling PsychPortAudio('Open', ...) will fail.

######################################################################
Stimulus timing:
Since each sound is played separately (not as a sequence of 5) care must be taken with the sound output latency. Our system has an approximate output latency of 10 ms (as indicated by psychtoolbox after calling PsychPortAudio('Open', ...) ). This is taken into account on line 163 where a pause of 140 ms is initiated. Together with the output latency it will result in approximately 150 ms SOA. The best way to check whether the SOAs are correct is to run the sequence for a while and then look at the output of diff(tz). It should be a repeating requence of aprox. 0.15 0.15 0.15 0.15 0.8-1.1. Some jitter of a few ms is unfortunately unavoidable on a standard computer.
 
######################################################################
Triggers:
I like to have a trigger at the onset of each sound. We have a dedicated trigger box that sends these triggers and the respective code is commented out on lines 161-162.
This code just has to be replaced by anything that works for your system. For example, the putvalue(parport, ...) lines for controlling the parallel port.
Don't forget to add any configuration code for the trigger device befor you start sending triggers.

######################################################################
Pauses:
At the moment there is a 5-second pause after each block (line 173). Only after the 4th block the KbStrokeWait function is called and a keypress (any key) is necessary to continue running the code (line 171). I use this opportunity to start a new EEG file, but pauses can be changed according to your own wishes because neither Bekinschtein nor Faugeras specify any pauses at all. 

######################################################################
Differences from what's reported in Faugeras et al. (2012):
1. There is an amplitude modulations for the 2. and 3. partials of a chord. So the normalized amplitudes for the three partials are 1, 1/2 and 1/4, respectively. It sounds more harmonic and humane. 
2. It is not reported in Faugeras et al. (2012), but I have made sure that there are at least 2 standards between any 2 deviants. Otherwise P3b may be attenuated. 



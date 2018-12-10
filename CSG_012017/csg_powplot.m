function csg_powplot(hand,handles)

% FORMAT csg_hypnoplot(hand,handles,windowsize,score,scorer)
% Plotting hypnogram for scored sleep data set.
%__________________________________________________________________
% Copyright (C) 2016 Coppieters 

% Written by Coppieters 2016
% Coma Science Group, University of Liège, Belgium
% $Id$

set(handles.figure1,'CurrentAxes',hand)
cla(hand);

totime = ceil(nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1}));
 xtick =sort((totime):-(totime/10):0);
%     set(handles.axes4,'Xtick',sort((totime):-(totime/10):0));
    if isfield(handles.Dmeg{1},'info') &&  isfield(handles.Dmeg{1}.info,'hour')
            xtick = mod(xtick + handles.offset,24*60^2);
    elseif isfield(handles,'offset')
        xtick = mod(xtick + handles.offset,24*60^2);
    end
    
    [dumb, times] = crc_time_converts(xtick);
    hold on,
    set(handles.axes4,'YTick',[0:10:floor(handles.pow.frequency(end))], ...
        'YTickLabel', [0:10:floor(handles.pow.frequency(end))],...
        'Xtick',xtick,...
        'XTickLabel',times, ...
        'Fontsize',8);
    surf(handles.pow.tempo+handles.powinfo.epoch/2 + xtick(1),handles.pow.frequency,handles.pow.power','EdgeColor','none')
    colormap(jet);
    xlim([xtick(1) xtick(end)])
ylim([0 floor(handles.pow.frequency(end))])
% 
set(handles.figure1,'CurrentAxes',handles.axes1)
if isfield(handles,'slider1')
    hold on
    slidval = get(handles.slider1,'Value');
    pos = min(slidval+handles.winsize/2,nsamples(handles.Dmeg{1})/fsample(handles.Dmeg{1}));
    handles.cursor = plot(pos,0,'^','Color',[0.2 0.2 0.2], ...
        'LineWidth',2.5,'tag','cursor');
    hold off
end

return

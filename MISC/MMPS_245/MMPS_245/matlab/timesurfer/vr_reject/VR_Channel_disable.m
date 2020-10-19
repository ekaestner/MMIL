function VR_Channel_disable(h,~)   %Allows for disabling of channels by clicking title objects
global TRIALDATA
global CHANNELDATA
if ~strcmp('All',CHANNELDATA.view_mode) && ~strcmp('Inspection',CHANNELDATA.view_mode)
%     chan_lab=length(CHANNELDATA.label{1});
    tmp_str                        = get(h,'String');
    tmp_string                     = strtok(tmp_str,'-');
    chan_index                     = find(strcmp(tmp_string,CHANNELDATA.label(:)));
    CHANNELDATA.toggle(chan_index) = 0;
    
    set(CHANNELDATA.patch(chan_index),'Visible','off','HitTest','off');
    
    for cc=1:length(CHANNELDATA.plotdata)
        set(CHANNELDATA.plotdata(cc).threshold_lines(chan_index,1),'Visible','off','HitTest','off');
        set(CHANNELDATA.plotdata(cc).threshold_lines(chan_index,2),'Visible','off','HitTest','off');
        tmp_title = get(CHANNELDATA.plotdata(cc).axes(chan_index),'Title');
        set(tmp_title,'Color',[0 0 1],'ButtonDownFcn',{@VR_Channel_enable});
        set(CHANNELDATA.plotdata(cc).plot_lines(chan_index,:),'Visible','off','HitTest','off');
    end;
end;
chan_data  = ~CHANNELDATA.toggle;
trial_data = ~TRIALDATA.toggle;
set(CHANNELDATA.badpanel,'Data',CHANNELDATA.label(chan_data))
set(TRIALDATA.badpanel,'Data',TRIALDATA.trial_cc(trial_data,:));
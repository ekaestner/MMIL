function VR_Channel_enable(h,~)
global TRIALDATA
global CHANNELDATA
if ~strcmp('All',CHANNELDATA.view_mode) && ~strcmp('Inspection',CHANNELDATA.view_mode)
    tmp_str    = get(h,'String');
    tmp_string = strtok(tmp_str,'-');
    chan_index=find(strcmp(tmp_string,CHANNELDATA.label(:)));
    set(CHANNELDATA.patch(chan_index),'Visible','on','HitTest','on');
    CHANNELDATA.toggle(chan_index)=1;
    for cc=1:length(CHANNELDATA.plotdata)
        set(CHANNELDATA.plotdata(cc).threshold_lines(chan_index,1),'Visible','on','HitTest','on');
        set(CHANNELDATA.plotdata(cc).threshold_lines(chan_index,2),'Visible','on','HitTest','on');
        tmp_title=get(CHANNELDATA.plotdata(cc).axes(chan_index),'Title');
        set(tmp_title,'Color',[0 0 0],'ButtonDownFcn',{@VR_Channel_disable});
        set(CHANNELDATA.plotdata(cc).plot_lines(chan_index,:),'Visible','on','HitTest','on');
    end;
end;
chan_data = ~CHANNELDATA.toggle;
trial_data = ~TRIALDATA.toggle;
set(CHANNELDATA.badpanel,'Data',CHANNELDATA.label(chan_data))
set(TRIALDATA.badpanel,'Data',TRIALDATA.trial_cc(trial_data,:));
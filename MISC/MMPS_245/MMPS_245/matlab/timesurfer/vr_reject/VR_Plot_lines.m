function VR_Plot_lines(h,~)
global CHANNELDATA
global TRIALDATA
[chan tmptrial] = find(CHANNELDATA.plotdata(CHANNELDATA.c_cond).plot_lines == h); %chan =channel index, tmptrial=condition trial number
trial = find(tmptrial == TRIALDATA.trial_cc(:,1) & CHANNELDATA.c_cond == TRIALDATA.trial_cc(:,2)); %trial= absolute trial number
max = TRIALDATA.max(2).data(trial,1);



if (strcmp('alt',get(gcf,'SelectionType'))) && ~(strcmp(CHANNELDATA.view_mode,'All')) && ~(strcmp(CHANNELDATA.view_mode,'Selection Mode')) && ~(strcmp(CHANNELDATA.view_mode,'Inspection'))    

    set(CHANNELDATA.plotdata(CHANNELDATA.c_cond).plot_lines(:,tmptrial),'Visible','off','HitTest','off');
    set(TRIALDATA.patch(trial),'Visible','off','HitTest','off');
    TRIALDATA.toggle(trial) = false;
    %code for updating threshold lines 
%     tmp           = find(TRIALDATA.toggle);
%     tmp_good_val  = TRIALDATA.max(2).data(tmp,1);
%     mean_all_chan = mean(tmp_good_val);
%     std_all_chan  = std(tmp_good_val);
%     set(TRIALDATA.oline,'YData',[mean_all_chan+2.5*std_all_chan mean_all_chan+2.5*std_all_chan]);
%     set(TRIALDATA.eoline,'YData',[mean_all_chan+5*std_all_chan mean_all_chan+5*std_all_chan]);

elseif(strcmp('Selection Mode',CHANNELDATA.view_mode))
 
    if (TRIALDATA.toggle(trial))
        
        lines = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).plot_lines,'Color',[1 0 0]);
        set(CHANNELDATA.plotdata(CHANNELDATA.c_cond).plot_lines(:,tmptrial),'Color',[1 0 0],'LineWidth',2);
        set(TRIALDATA.patch(trial),'Selected','on','HitTest','on');
        set(lines,'Color',[0 0 1]);
        TRIALDATA.toggle(trial) = false;
        
        
        
    
    elseif(~TRIALDATA.toggle(trial))
    
        set(CHANNELDATA.plotdata(CHANNELDATA.c_cond).plot_lines(:,tmptrial),'Color',[0 0 0],'LineWidth',.5);
        set(TRIALDATA.patch(trial),'Selected','off','HitTest','off');
        TRIALDATA.toggle(trial) = true;
  
    end;
    %code for updating threshold lines 
%     tmp           = find(TRIALDATA.toggle);
%     tmp_good_val  = TRIALDATA.max(2).data(tmp,1);
%     mean_all_chan = mean(tmp_good_val);
%     std_all_chan  = std(tmp_good_val);
%     set(TRIALDATA.oline,'YData',[mean_all_chan+2.5*std_all_chan mean_all_chan+2.5*std_all_chan]);
%     set(TRIALDATA.eoline,'YData',[mean_all_chan+5*std_all_chan mean_all_chan+5*std_all_chan]);
    
elseif (strcmp('Inspection',CHANNELDATA.view_mode))
    VR_Sum_points(TRIALDATA.patch(trial),[]); 
else
    tmp_lines     = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).plot_lines,'flat','Color',[0 0 1]);
    tmp_patch_c   = findobj(CHANNELDATA.patch,'flat','Selected','on');
    tmp_patch_t   = findobj(TRIALDATA.patch,'flat','Selected','on');
    tmp_raw_lines = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_lines,'flat','Visible','on');
    
    set(tmp_patch_c,'Selected','off');
    set(tmp_patch_t,'Selected','off');
    set(tmp_lines,'LineWidth',.05,'Color',[0 0 0],'ButtonDownFcn',{@VR_Plot_lines});
    set(tmp_raw_lines,'Visible','off');
    set(CHANNELDATA.plotdata(CHANNELDATA.c_cond).plot_lines(:,tmptrial),'LineWidth',2,'Color',[0 0 1],'ButtonDownFcn',{@VR_Hide_line});
    set(TRIALDATA.patch(trial),'Selected','on');
    set(CHANNELDATA.patch(chan),'Selected','on');
    set(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_lines(:,tmptrial),'Visible','on');
    set(CHANNELDATA.raw_plot,'Name',['Raw Plot showing Trial: ',num2str(tmptrial), ' in Condition Number: ',num2str(CHANNELDATA.c_cond)]);
    
    chan    = TRIALDATA.max(2).data(trial,2);
%     dis_str = ['Trial #: ' num2str(tmptrial),' Cond: ',num2str(CHANNELDATA.c_cond),' Max: ', num2str(max), ' on Channel Number: ' CHANNELDATA.label{chan}];
%     dis_str = [' Cond: ',num2str(CHANNELDATA.c_cond),'Trial #: ' num2str(trial),'- Max: ', num2str(max), ' on Channel Number: ' CHANNELDATA.label{chan}];
    dis_str = [' Cond ',num2str(CHANNELDATA.c_cond),', Trial ' num2str(trial),'- Max Amp ', num2str(max), ' on Channel ' CHANNELDATA.label{chan}];
    set(CHANNELDATA.print_txt,'String',dis_str);
    
    VR_Sum_points(CHANNELDATA.patch(chan),[],0);
end;
set(CHANNELDATA.badpanel,'Data',CHANNELDATA.label(~(CHANNELDATA.toggle)));
set(TRIALDATA.badpanel,'Data',TRIALDATA.trial_cc(~(TRIALDATA.toggle),:));

function VR_Hide_line(h,~)
global CHANNELDATA
global TRIALDATA
[chan tmptrial] = find(CHANNELDATA.plotdata(CHANNELDATA.c_cond).plot_lines == h);
trial=find(tmptrial == TRIALDATA.trial_cc(:,1) & CHANNELDATA.c_cond == TRIALDATA.trial_cc(:,2)); 
set(CHANNELDATA.plotdata(CHANNELDATA.c_cond).plot_lines(:,tmptrial),'LineWidth',.05,'Color',[0 0 0],'ButtonDownFcn',{@VR_Plot_lines});
set(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_lines(:,tmptrial),'Visible','off');
set(CHANNELDATA.patch(chan),'Selected','off')
chan_data = find(~CHANNELDATA.toggle);
trial_data = find(~TRIALDATA.toggle);
set(CHANNELDATA.badpanel,'Data',CHANNELDATA.label(chan_data))
set(TRIALDATA.badpanel,'Data',TRIALDATA.trial_cc(trial_data,:));
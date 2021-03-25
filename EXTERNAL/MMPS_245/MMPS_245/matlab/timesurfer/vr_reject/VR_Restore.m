function VR_Restore(h,index,type)
global TRIALDATA
global CHANNELDATA
data=get(h,'Data');

if type==0
    if ~isempty(index.Indices)
        % Restores Channel
        id = data{index.Indices(1)}; %Channel Label
        index = find(strcmp(id,CHANNELDATA.label)); %Channel Numerical Index
        tmp_string=get(CHANNELDATA.plotdata(1).axes(index),'Title');
        
        if(strcmp(CHANNELDATA.view_mode,'Inspection'))
            chan = CHANNELDATA.patch(index);
            VR_Sum_points(chan,[]);
        else
            VR_Channel_enable(tmp_string,[]);
            %updates threshold lines on summary plot after restoring channel
            tmp = CHANNELDATA.toggle(:);
            tmp_good_val = CHANNELDATA.max(2).data(tmp,1);
            mean_all_chan = mean(tmp_good_val);
            std_all_chan=std(tmp_good_val);
            set(CHANNELDATA.oline,'XData',[mean_all_chan+2.5*std_all_chan mean_all_chan+2.5*std_all_chan]);
            set(CHANNELDATA.eoline,'XData',[mean_all_chan+5*std_all_chan mean_all_chan+5*std_all_chan]);
        end;
    end;
elseif type==1
    if ~isempty(index.Indices)
        if ~iscell(data); data=num2cell(data); end;
        id_trial = data{index.Indices(1),1};
        id_cond = data{index.Indices(1),2};
        %         id_trial=str2num(id_trial);
        %         id_cond=str2num(id_cond);
        if ischar(id_trial); id_trial=str2num(id_trial); end;
        if ischar(id_cond); id_cond=str2num(id_cond); end;
        indx = find((id_trial==TRIALDATA.trial_cc(:,1)) & (id_cond==TRIALDATA.trial_cc(:,2))); %absolute trial num
        if (strcmp(CHANNELDATA.view_mode,'Inspection'))
            trial = TRIALDATA.patch(indx);
            VR_Sum_points(trial,[]);
        else
            TRIALDATA.toggle(indx)=true;
            set(TRIALDATA.patch(indx),'Visible','on','Selected','off','HitTest','on','ButtonDownFcn',{@VR_Sum_points});
            list=find(CHANNELDATA.toggle);
            tmp=CHANNELDATA.plotdata(id_cond).plot_lines(list,id_trial);
            set(tmp,'Visible','on','Color',[0 0 0],'Selected','off','HitTest','on','ButtonDownFcn',{@VR_Plot_lines},'LineWidth',.5);
            % Update Trial Threshold lines following restore
            tmp= find(TRIALDATA.toggle);
            tmp_good_val = TRIALDATA.max(2).data(tmp,1);
            mean_all_chan = mean(tmp_good_val);
            std_all_chan=std(tmp_good_val);
            set(TRIALDATA.oline,'YData',[mean_all_chan+2.5*std_all_chan mean_all_chan+2.5*std_all_chan])
            set(TRIALDATA.eoline,'YData',[mean_all_chan+5*std_all_chan mean_all_chan+5*std_all_chan])
        end;
    end;
end;
chan_data = CHANNELDATA.label(~CHANNELDATA.toggle);
trial_data = ~TRIALDATA.toggle;
temp = TRIALDATA.trial_cc(trial_data,:);
% T_data = {};
if ~isempty(temp)
    T_data(:,1) = num2cell(temp(:,1));
    T_data(:,2) = num2cell(temp(:,2));
else
    T_data={};
end;
set(TRIALDATA.badpanel,'Data',T_data);
set(CHANNELDATA.badpanel,'Data',chan_data)
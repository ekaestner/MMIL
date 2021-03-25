function VR_Sum_points(varargin)%(h,x,code)
h=varargin{1};
if length(varargin)>=3
    code=varargin{3};
    %disp('run 2');
else
    code=1; %first pass
    %disp('run 1');
end;

global CHANNELDATA
global TRIALDATA
tmp_val=get(h,'tag');
i = 0;
if strcmp(tmp_val,'chan')
    type=0; %Channel
    i = find(CHANNELDATA.patch(:)==h); %Chan index
    name=CHANNELDATA.label(i); %Channel Name
%     if(strcmp(CHANNELDATA.view_mode,'All'))
    if(strcmp(CHANNELDATA.view_mode,'All')) || (strcmp(CHANNELDATA.view_mode,'Selection Mode')) || (strcmp(CHANNELDATA.view_mode,'Inspection'))
        value=CHANNELDATA.max(1).data(i,1); %Max Value
    else
        value=CHANNELDATA.max(2).data(i,1); %Max Value
    end;
else
    type=1; %Trial
    i = find(TRIALDATA.patch(:)==h); %abs trial index
    name = TRIALDATA.trial_cc(i,1); %cond trial index
    cc = TRIALDATA.trial_cc(i,2);
    if(strcmp(CHANNELDATA.view_mode,'All')) || (strcmp(CHANNELDATA.view_mode,'Selection Mode')) || (strcmp(CHANNELDATA.view_mode,'Inspection'))
        value = TRIALDATA.max(1).data(i,1);
        chan = TRIALDATA.max(1).data(i,2);
    else
        value = TRIALDATA.max(2).data(i,1);
        chan = TRIALDATA.max(2).data(i,2);
    end;
    if code==3
        chan=varargin{4};
    end;
end;

if (strcmp('alt',get(gcf,'SelectionType'))) && (~strcmp(CHANNELDATA.view_mode,'All')) && (~strcmp(CHANNELDATA.view_mode,'Inspection')) || (strcmp(CHANNELDATA.view_mode,'Selection Mode'));
    %%remove points
    if type == 0  %channels
        if ~strcmp(CHANNELDATA.view_mode,'Selection Mode') %Turn Patch off outside of selection mode
            set(h,'Visible','off','HitTest','off');
        else
            set(h,'HitTest','on','Selected','on','ButtonDownFcn',{@VR_hidepoint});
        end;
        CHANNELDATA.toggle(i) = 0;
%         % Does not update outlier lines in selection mode.
%         if ~(strcmp(CHANNELDATA.view_mode,'Selection Mode'))
%             tmp = CHANNELDATA.toggle(:);
%             tmp_good_val = CHANNELDATA.max(2).data(tmp,1);
%             mean_all_chan = mean(tmp_good_val);
%             std_all_chan=std(tmp_good_val);
%             set(CHANNELDATA.oline,'XData',[mean_all_chan+2.5*std_all_chan mean_all_chan+2.5*std_all_chan]);
%             set(CHANNELDATA.eoline,'XData',[mean_all_chan+5*std_all_chan mean_all_chan+5*std_all_chan]);
%         end;
%         
        for cc1=1:length(CHANNELDATA.plotdata)
            set(CHANNELDATA.plotdata(cc1).threshold_lines(i,:),'Visible','off','HitTest','off');
            tmp_title=get(CHANNELDATA.plotdata(cc1).axes(i),'Title');
            set(tmp_title,'Color',[0 0 1]);
            tmp=CHANNELDATA.plotdata(cc1).plot_lines(i,:);
            set(tmp,'Visible','off','HitTest','off');
        end;
    else
        
        TRIALDATA.toggle(i)=0;
        if ~strcmp(CHANNELDATA.view_mode,'Selection Mode')
            set(h,'Visible','off','HitTest','off');
        else
            set(h,'HitTest','on','Selected','on','ButtonDownFcn',{@VR_hidepoint});
        end;
%         if ~(strcmp(CHANNELDATA.view_mode,'Selection Mode'))
%             tmp= find(TRIALDATA.toggle);
%             tmp_good_val = TRIALDATA.max(2).data(tmp,1);
%             mean_all_chan = mean(tmp_good_val);
%             std_all_chan=std(tmp_good_val);
%             set(TRIALDATA.oline,'YData',[mean_all_chan+2.5*std_all_chan mean_all_chan+2.5*std_all_chan])
%             set(TRIALDATA.eoline,'YData',[mean_all_chan+5*std_all_chan mean_all_chan+5*std_all_chan])
%         end;
        
        %remove lines
        if ~strcmp(CHANNELDATA.view_mode,'Selection Mode')
            set(CHANNELDATA.plotdata(cc).plot_lines(:,name),'Visible','off','HitTest','on')
        else
            lines = findobj(CHANNELDATA.plotdata(cc).plot_lines,'Color',[1 0 0]);
            set(CHANNELDATA.plotdata(cc).plot_lines(:,name),'Color',[1 0 0],'ButtonDownFcn',{@VR_Plot_lines},'LineWidth',2);
            set(lines,'Color',[0 0 1]);
        end;
    end;
    % elseif (strcmp(CHANNELDATA.view_mode,'Inspection'))
    %     if type == 0
    %         dis_str = ['Channel: ',name{:},' Max: ', num2str(value), ' on Trial: ' num2str(CHANNELDATA.max(2).data(i,2)),' of condition ' num2str(CHANNELDATA.max(2).data(i,3))];
    %         if code==1 || code == 3
    %         set(CHANNELDATA.print_txt,'String',dis_str);
    %         end;
    %         set(h,'Selected','on');
    %         set(h,'ButtonDownFcn',{@VR_hidepoint});
    %         %Select trials
    %         abst = find(TRIALDATA.trial_cc(:,1)==CHANNELDATA.max(2).data(i,2) & TRIALDATA.trial_cc(:,2)==CHANNELDATA.max(2).data(i,3));
    %         if (strcmp(get(TRIALDATA.patch(abst),'Visible'),'on') && code==1)
    %             %disp('Hello World the next step is wrong');
    %             VR_Sum_points(TRIALDATA.patch(abst),[],3);
    %         end;
    %         %change page
    %         if code == 1 || code == 3
    %         pan_h = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels,'flat','Visible','on');
    %         set(pan_h,'Visible','off');
    %         new_h = get(CHANNELDATA.plotdata(CHANNELDATA.c_cond).axes(i),'Parent');
    %         set(new_h,'Visible','on');
    %         end;
    %     elseif type==1
    %
    %     end;
    %     %Insert inspection mode code
    %
else
    %% show points
    if type == 0
        tmpobj = findobj(CHANNELDATA.patch(:),'Selected','on');
        set(tmpobj,'Selected','off','ButtonDownFcn',{@VR_Sum_points});
        tmpobj = findobj(CHANNELDATA.patch(:),'EdgeColor',[1 0 0]);
        if (~strcmp(CHANNELDATA.view_mode,'Inspection'))
            set(tmpobj,'EdgeColor',[0 0 0]);
        else
            set(tmpobj,'EdgeColor',[0 0 0]);
            
        end;
        %         fprintf('Channel: %s  Max: %3.2e \n',name{:},value);
        if (strcmp(CHANNELDATA.view_mode,'All')) || (strcmp(CHANNELDATA.view_mode,'Selection Mode')) || (strcmp(CHANNELDATA.view_mode,'Inspection'))
            dis_str = ['Channel: ',name{:},' Max: ', num2str(value), ' on Trial: ' num2str(CHANNELDATA.max(1).data(i,2)),' of condition ' num2str(CHANNELDATA.max(1).data(i,3))];
        else
            dis_str = ['Channel: ',name{:},' Max: ', num2str(value), ' on Trial: ' num2str(CHANNELDATA.max(2).data(i,2)),' of condition ' num2str(CHANNELDATA.max(2).data(i,3))];
        end;
        if code==1
            set(CHANNELDATA.print_txt,'String',dis_str);
            
        end;
        
        if (~strcmp(CHANNELDATA.view_mode,'Inspection'))
            set(h,'Selected','on');
            set(h,'ButtonDownFcn',{@VR_hidepoint});
        else
            tmp=findobj(TRIALDATA.patch,'Selected','on');
            set(tmp,'Selected','off');
            set(h,'Selected','on','EdgeColor',[1 0 0]);
        end;
        
        %change page
        if code == 1
            pan_h = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels,'flat','Visible','on');
            set(pan_h,'Visible','off');
            new_h = get(CHANNELDATA.plotdata(CHANNELDATA.c_cond).axes(i),'Parent');
            set(new_h,'Visible','on');
            raw = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan,'flat','Visible','on');
            set(raw,'Visible','off');
            raw_new = get(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_axes(i),'Parent');
            set(raw_new,'Visible','on');
        end;
        %Select trials
        if strcmp('Inspection',CHANNELDATA.view_mode)
            abst = find(TRIALDATA.trial_cc(:,1)==CHANNELDATA.max(1).data(i,2) & TRIALDATA.trial_cc(:,2)==CHANNELDATA.max(1).data(i,3));          
        else
            abst = find(TRIALDATA.trial_cc(:,1)==CHANNELDATA.max(2).data(i,2) & TRIALDATA.trial_cc(:,2)==CHANNELDATA.max(2).data(i,3));
        end;
%         abst = find(TRIALDATA.trial_cc(:,1)==CHANNELDATA.max(2).data(i,2) & TRIALDATA.trial_cc(:,2)==CHANNELDATA.max(2).data(i,3));
        if (strcmp(get(TRIALDATA.patch(abst),'Visible'),'on') && code==1)
            %disp('Hello World the next step is wrong');
            if strcmp('Inspection',CHANNELDATA.view_mode)
                VR_Sum_points(TRIALDATA.patch(abst),[],3,i);
            else
                VR_Sum_points(TRIALDATA.patch(abst),[],0);
            end;
        end;
        
        %Jump  conditions
    end;
    if type == 1
        if (~strcmp(CHANNELDATA.view_mode,'Inspection'))
            tmpobj = findobj(TRIALDATA.patch(:),'Selected','on');
            set(tmpobj,'Selected','off','ButtonDownFcn',{@VR_Sum_points});
        end;
        x = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).plot_lines,'flat','Color',[0 0 1]); %finding all the blue lines(selected)
        %Change raw lines, change line colors for appropriate mode
        for cci=1:CHANNELDATA.num_conds
            tmp_raw = findobj(CHANNELDATA.plotdata(cci).raw_lines,'flat','Visible','on');
            set(tmp_raw,'Visible','off');
            if strcmp(CHANNELDATA.view_mode,'Inspection')
                x = findobj(CHANNELDATA.plotdata(cci).plot_lines,'flat','Color',[1 0 0]);
                set(x,'Color',[0 0 1]);
            else
                x = findobj(CHANNELDATA.plotdata(cci).plot_lines,'flat','Color',[0 0 1]);
                set(x,'LineWidth',.05,'Color',[0 0 0]);
            end;
        end;
        
        if (strcmp(get(CHANNELDATA.patch(chan),'Visible'),'on') && code==1)
            %disp('Wrong');
            if (strcmp('Inspection',CHANNELDATA.view_mode))
                VR_Sum_points(CHANNELDATA.patch(chan),[],3);
            else
                VR_Sum_points(CHANNELDATA.patch(chan),[],0);
            end;
            
        end;
        %         if (~strcmp(CHANNELDATA.view_mode,'Inspection'))
        %             tmpobj = findobj(CHANNELDATA.patch(:),'Selected','on');
        %             set(tmpobj,'Selected','off','ButtonDownFcn',{@VR_Sum_points});
        %         end;
        
        
        tmpobj = findobj(TRIALDATA.patch(:),'EdgeColor',[1 0 0]);
        set(tmpobj,'EdgeColor',[0 0 0]);
        
        
        if code==1
%             dis_str = [' Cond: ',num2str(cc),'Trial #: ' num2str(name),'- Max: ', num2str(value), ' on Channel Number: ' CHANNELDATA.label{chan}];
            dis_str = [' Cond ',num2str(cc),', Trial ' num2str(name),'- Max Amp ', num2str(value), ' on Channel ' CHANNELDATA.label{chan}];
            set(CHANNELDATA.print_txt,'String',dis_str);
        end;
        
        if ~(strcmp(CHANNELDATA.view_mode,'Inspection'))
            set(h,'Selected','on');
            set(h,'ButtonDownFcn',{@VR_hidepoint});
            set(CHANNELDATA.plotdata(cc).plot_lines(:,name),'LineWidth',2,'Color',[0 0 1]);
        else
            set(h,'Selected','on','EdgeColor',[1 0 0]);
            set(CHANNELDATA.plotdata(cc).plot_lines(:,name),'LineWidth',2,'Color',[1 0 0]);
        end;
        
        set(CHANNELDATA.plotdata(cc).raw_lines(:,name),'Visible','on');
        
        if strcmp(CHANNELDATA.view_mode,'Inspection')
            cc_h=findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels,'flat','Visible','on');
            cc_h_raw = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan,'flat','Visible','on');
            set(CHANNELDATA.text(CHANNELDATA.c_cond),'Visible','off');
            CHANNELDATA.c_cond = cc;
            set(CHANNELDATA.text(CHANNELDATA.c_cond),'Visible','on');
            new_h = get(CHANNELDATA.plotdata(CHANNELDATA.c_cond).axes(chan),'Parent');
            new_h_raw = get(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_axes(chan),'Parent');
            set(cc_h,'Visible','off');
            axes = get(cc_h,'Children');
%             set(axes,'Ylim',CHANNELDATA.Ylim);
            set(new_h,'Visible','on');
            set(cc_h_raw,'Visible','off');
            axes_raw = get(cc_h_raw,'Children');
            set(axes_raw,'Ylim',CHANNELDATA.Ylim_raw);
            set(new_h_raw,'Visible','on');
            set(CHANNELDATA.plotfig,'Name',['Condition: ',num2str(CHANNELDATA.c_cond), ' View Mode: ',CHANNELDATA.view_mode]);
            tmp_raw_lines = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_lines,'Visible','on');
            if ~isempty(tmp_raw_lines)
                [raw_info(1) raw_info(2)] = find(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_lines == tmp_raw_lines(1));
                set(CHANNELDATA.raw_plot,'Name',['Raw Plot - in Condition # ',num2str(CHANNELDATA.c_cond), ' Trial # ',num2str(raw_info(2))]);
            else
                set(CHANNELDATA.raw_plot,'Name',['Raw Plot']);
            end;
        end;
        
        if code==1 || code == 3
            set(CHANNELDATA.raw_plot,'Name',['Raw Plot showing ', num2str(name),'on Condition ', num2str(cc)]);
        end;
        
        if (strcmp('alt',get(gcf,'SelectionType'))) && code == 1 && strcmp(CHANNELDATA.view_mode,'Inspection')
            if type==1
                if TRIALDATA.toggle(i)
                    TRIALDATA.toggle(i)       = 0;
                else
                    TRIALDATA.toggle(i)       = 1;
                end;
            elseif type==0;
                if CHANNELDATA.toggle(i)
                    CHANNELDATA.toggle(i) = 0;
                else
                    CHANNELDATA.toggle(i) = 1;
                end;
            end;
        end;
        
        
        
    end;
end;
chan_data = CHANNELDATA.label(~CHANNELDATA.toggle);
trial_data = find(~TRIALDATA.toggle);
temp = TRIALDATA.trial_cc(trial_data,:);
T_data(:,1) = cellstr(num2str(temp(:,1)));
T_data(:,2) = cellstr(num2str(temp(:,2)));
%Remove old(the list is reset everytime so no need to remove old) set current if code==1
set(TRIALDATA.badpanel,'Data',T_data);
set(CHANNELDATA.badpanel,'Data',chan_data)
clear h x temp data ind tmp_str_i

function VR_hidepoint(h,x)
global CHANNELDATA
global TRIALDATA
tmp_val=get(h,'tag');
i = 0;
%Determine if patch selected was a channel or a trial patch
if strcmp(tmp_val,'chan')
    type=0; %Channel
    i = find(CHANNELDATA.patch(:)==h);
    name=CHANNELDATA.label(i);
    value=CHANNELDATA.max(2).data(i,1);
else
    type=1; %Trial
    i = find(TRIALDATA.patch(:)==h);
    name = TRIALDATA.trial_cc(i,1);
    cc = TRIALDATA.trial_cc(i,2);
    value = TRIALDATA.max(2).data(i,1);
    chan = TRIALDATA.max(2).data(i,2);
end;


if (strcmp(CHANNELDATA.view_mode,'Selection Mode'))
    set(h,'HitTest','on','Selected','off','ButtonDownFcn',{@VR_Sum_points});
    if type == 0 %Channel
        CHANNELDATA.toggle(i) = logical(1);
        for cc=1:length(CHANNELDATA.plotdata)
            set(CHANNELDATA.plotdata(cc).threshold_lines(i,:),'Visible','on','HitTest','off');
            tmp_title=get(CHANNELDATA.plotdata(cc).axes(i),'Title');
            set(tmp_title,'Color',[0 0 0]);
            tmp=CHANNELDATA.plotdata(cc).plot_lines(i,:);
            set(tmp,'Visible','on','HitTest','on');
        end;
    elseif type == 1 %Trial
        TRIALDATA.toggle(i) = logical(1);
        set(CHANNELDATA.plotdata(cc).plot_lines(:,name),'Color',[0 0 0],'ButtonDownFcn',{@VR_Plot_lines});
        pan_h = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels,'flat','Visible','on');
        set(pan_h,'Visible','off');
        new_h = get(CHANNELDATA.plotdata(CHANNELDATA.c_cond).axes(i),'Parent');
        set(new_h,'Visible','on');
    end;
else
    if type == 0
        dis_str = ['Channel: ',name{:},' Max: ', num2str(value), ' on Trial: ' num2str(CHANNELDATA.max(2).data(i,2)),' of condition ' num2str(CHANNELDATA.max(2).data(i,3))];
        set(CHANNELDATA.print_txt,'String',dis_str);
        set(h,'Selected','off')
        set(h,'ButtonDownFcn',{@VR_Sum_points});
    end;
    if type == 1
        dis_str = [' Cond ',num2str(cc),', Trial ' num2str(name),'- Max Amp ', num2str(value), ' on Channel ' CHANNELDATA.label{chan}];
        set(CHANNELDATA.print_txt,'String',dis_str);
        set(h,'Selected','off');
        set(h,'ButtonDownFcn',{@VR_Sum_points});
        x=findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).plot_lines,'flat','Color',[0 0 1]);
        set(x,'LineWidth',.05,'Color',[0 0 0]);
    end;
end;
chan_data = CHANNELDATA.label(~CHANNELDATA.toggle);
trial_data = find(~TRIALDATA.toggle);
temp = TRIALDATA.trial_cc(trial_data,:);
data(:,1) = cellstr(num2str(temp(:,1)));
data(:,2) = cellstr(num2str(temp(:,2)));
set(TRIALDATA.badpanel,'Data',data);
set(CHANNELDATA.badpanel,'Data',chan_data)
clear h x
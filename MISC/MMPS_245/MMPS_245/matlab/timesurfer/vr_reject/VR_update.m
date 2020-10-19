% function VR_update(varargin)
% varargin
function VR_update(h,x,epoch_data,tidx)
global CHANNELDATA
global TRIALDATA
chan_data = find(~CHANNELDATA.toggle);
chan_good = find(CHANNELDATA.toggle);
trial_data = find(~TRIALDATA.toggle);
cc_bad=TRIALDATA.trial_cc(trial_data,2);
trial_cc_bad=TRIALDATA.trial_cc(trial_data,1);
trial_good=find(TRIALDATA.toggle);



for ci=1:epoch_data.num_sensors
    for cc=1:CHANNELDATA.num_conds
        tmp_trial=find(TRIALDATA.toggle & TRIALDATA.trial_cc(:,2)==1);
        list=TRIALDATA.trial_cc(tmp_trial,1);
        if cc == 1
            [CHANNELDATA.max(2).data(ci,1) CHANNELDATA.max(2).data(ci,2)] = max(max(abs(epoch_data.epochs(1).data(ci,tidx,list)),[],2));
            CHANNELDATA.max(2).data(ci,3) = 1;
            h=CHANNELDATA.patch(ci); %set xdata
            set(h,'XData',CHANNELDATA.max(2).data(ci,1));
        else
            tmp_trial=find(TRIALDATA.toggle & TRIALDATA.trial_cc(:,2)==cc);
            list=TRIALDATA.trial_cc(tmp_trial,1);
            [x y] = max(max(abs(epoch_data.epochs(cc).data(ci,tidx,list)),[],2));
            y=list(y);
            if x > CHANNELDATA.max(2).data(ci,1)
                CHANNELDATA.max(2).data(ci,1) = x;
                CHANNELDATA.max(2).data(ci,2) = y;
                CHANNELDATA.max(2).data(ci,3) = cc;
                h=CHANNELDATA.patch(ci);
                set(h,'XData',CHANNELDATA.max(2).data(ci,1));
            end;
        end;
    end
end;



ind=1;
for cc=1:length(epoch_data.epochs)
    for ti=1:epoch_data.epochs(cc).num_trials
        [TRIALDATA.max(2).data(ind,1) y] = max(max(abs(epoch_data.epochs(cc).data(chan_good,tidx,ti)),[],2)); %max val, channel
        TRIALDATA.max(2).data(ind,2) = chan_good(y);
        set(TRIALDATA.patch(ind),'Ydata',TRIALDATA.max(2).data(ind,1));
        ind=ind+1;
    end;
end;
clear ind

%----
%determine bad stuff
ncond  = length(CHANNELDATA.plotdata);
nchan  = size(  CHANNELDATA.plotdata(cc).plot_lines,1);
for cc = 1:ncond
  ntrial = size(CHANNELDATA.plotdata(cc).plot_lines,2);
  
    for ci = 1:nchan
            thresdata1  = get(CHANNELDATA.plotdata(cc).threshold_lines(ci,1),'Ydata');
            m_thres1    = abs(mean(thresdata1));
            thresdata2  = get(CHANNELDATA.plotdata(cc).threshold_lines(ci,2),'Ydata');
            m_thres2    = abs(mean(thresdata2));
        for ti = 1:ntrial
            y_data = get(CHANNELDATA.plotdata(cc).plot_lines(ci,ti),'Ydata');

            trial_check = find((TRIALDATA.trial_cc(:,1)==ti) & (TRIALDATA.trial_cc(:,2)==cc));
            if (TRIALDATA.toggle(trial_check)) && (any(y_data > thresdata1) || any(y_data < thresdata2))
                CHANNELDATA.plotdata(cc).toggle(ci,ti)=logical(0);
            else
                CHANNELDATA.plotdata(cc).toggle(ci,ti)=logical(1);
            end;
        end;
        bad_list=find(~CHANNELDATA.plotdata(cc).toggle(ci,:));

        if ~isempty(bad_list)
        dist =[];
        for ti=1:length(bad_list)
            y_data = get(CHANNELDATA.plotdata(cc).plot_lines(ci,bad_list(ti)),'Ydata');
            [max_trial pos]=max(abs(y_data));
            if any(y_data > thresdata1)
                dist(end+1) = abs((y_data(pos)-thresdata1(pos)));
            elseif any(y_data < thresdata2)
                dist(end+1) = abs((y_data(pos)-thresdata2(pos)));
            end;
        end;
        [Y,i]=sort(dist,'descend');
        if numel(i)>3 % added to prevent crowding 6.6.11 -BQR
            i = i(1:3);
        end
        tmp_title=get(CHANNELDATA.plotdata(cc).axes(ci),'Title');
        set(tmp_title,'String',strcat(CHANNELDATA.label(ci),'-',num2str(bad_list(i))));
        set(CHANNELDATA.plotdata(cc).axes(ci),'Title',tmp_title);
        else     
            tmp_title=get(CHANNELDATA.plotdata(cc).axes(ci),'Title');
        set(tmp_title,'String',strcat(CHANNELDATA.label(ci)));
        end;
    end;
end;

num_total = [epoch_data.epochs.num_trials];
num_bad   = arrayfun(@(x) find(TRIALDATA.trial_cc(~TRIALDATA.toggle,2)==x),1:length(num_total),'UniformOutput',0);
num_bad   = cell2mat(arrayfun(@(x) length(num_bad{x}),1:length(num_total),'UniformOutput',0));
num_good = num_total-num_bad;


%re calc outlier lines
%trial
tmp_good_val = TRIALDATA.max(2).data(trial_good,1);
mean_all_trial = mean(tmp_good_val);
std_all_trial = std(tmp_good_val);
set(TRIALDATA.oline,'YData',[mean_all_trial+2.5*std_all_trial mean_all_trial+2.5*std_all_trial]);
set(TRIALDATA.eoline,'YData',[mean_all_trial+5*std_all_trial mean_all_trial+5*std_all_trial]);
%channel
tmp_good_val = CHANNELDATA.max(2).data(chan_good,1);
mean_all_chan = mean(tmp_good_val);
std_all_chan = std(tmp_good_val);
set(CHANNELDATA.oline,'XData',[mean_all_chan+2.5*std_all_chan mean_all_chan+2.5*std_all_chan]);
set(CHANNELDATA.eoline,'XData',[mean_all_chan+5*std_all_chan mean_all_chan+5*std_all_chan]);
%get values for outlier and extreme outliers
e_value = get(CHANNELDATA.eoline,'Xdata');
o_value = get(CHANNELDATA.oline,'Xdata');
for ci = 1:length(CHANNELDATA.patch)
    tmp = get(CHANNELDATA.patch(ci),'Xdata');
    if(CHANNELDATA.toggle(ci)) && any(tmp>e_value)
        CHANNELDATA.outliers(ci,1)=logical(1);
    else
        CHANNELDATA.outliers(ci,1)=logical(0);
    end;
    if any(tmp>o_value) && (CHANNELDATA.toggle(ci))
        CHANNELDATA.outliers(ci,2)=logical(1);
    else
        CHANNELDATA.outliers(ci,2)=logical(0);
    end;
end;

for ti=1:length(TRIALDATA.patch)
    e_value=get(TRIALDATA.eoline,'Ydata');
    o_value=get(TRIALDATA.oline,'Ydata');
    tmp=get(TRIALDATA.patch(ti),'Ydata');
    if any(tmp>e_value) && (TRIALDATA.toggle(ti))
        TRIALDATA.outliers(ti,1)=logical(1);
    else
        TRIALDATA.outliers(ti,1)=logical(0);
    end;
    if any(tmp>o_value) && (TRIALDATA.toggle(ti))
        TRIALDATA.outliers(ti,2)=logical(1);
    else
        TRIALDATA.outliers(ti,2)=logical(0);
    end;
end;
tot_text  = ['Total Trials: ' num2str(num_total)];
good_text = ['Good Trials : ' num2str(num_good)];
etext=['Extreme Outlier Threshold: ' num2str(e_value(1)) ' (*)'];
otext=['Outlier Threshold Value: ' num2str(o_value(1))];
disp('Summary of Outlier Data');
disp(tot_text);
disp(good_text);
disp(etext);
disp(otext);
set(CHANNELDATA.etext,'String',etext);
set(CHANNELDATA.otext,'String',otext);
for cc = 1:CHANNELDATA.num_conds
    textline = ['Condition ' num2str(cc) '(Event Code: ' num2str(epoch_data.epochs(cc).event_code) '): '];
    tmp = find(TRIALDATA.outliers(:,1));
    tmp_trial = find((TRIALDATA.trial_cc(tmp,2))==cc);
    e_list = TRIALDATA.trial_cc(tmp(tmp_trial),1);
    tmp = find(TRIALDATA.outliers(:,2));
    tmp_trial = find((TRIALDATA.trial_cc(tmp,2))==cc);
    o_list = TRIALDATA.trial_cc(tmp(tmp_trial),1);
    for ei = 1:length(e_list)
        textline = [textline, num2str(e_list(ei)),'* '];
    end;
    for oi = 1:length(o_list)
        if ~ismember(o_list(oi),e_list)
            textline=[textline, ' ' num2str(o_list(oi))];
        end;
    end;
    set(CHANNELDATA.text(cc),'String',textline);
    %drawnow;
    disp(textline);
end;
set(CHANNELDATA.text(CHANNELDATA.c_cond),'Visible','on');

    
    
chan_data = find(~CHANNELDATA.toggle);
trial_data = find(~TRIALDATA.toggle);
set(CHANNELDATA.badpanel,'Data',CHANNELDATA.label(chan_data))
set(TRIALDATA.badpanel,'Data',TRIALDATA.trial_cc(trial_data,:));
set(CHANNELDATA.print_txt,'String','');
disp('done');
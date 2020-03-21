%% Function to read events & create file to epoch NYU iEEG data
% Function to epoch NYU iEEG data. 
% 
% fcfg.dataset = path to continuous data structure you will be epoching
% fcfg.dsfact = rate that the samples must be divided by (if applicable) - EJK need to add
% fcfg.ignore = vector of trial indicies to ignore
% fcfg.minduration = smallest number of time between evenst (in seconds) 
% fcfg.timelim = 2 vectors of the time the events should be between (in seconds)
% fcfg.evt = vector of events to include in the final list
% fcfg.pre = time (in s) to include before the event
% fcfg.pos = time (in s) to include after the event
% fcfg.trg_chn = language to use for trg_chn
% 
% Created by Erik Kaestner (4-4-14) ekaestne@ucsd.edu

function [trl, evt] = ft_pedot_trialfun(fcfg)

%% Set up to find the timings
Fs = fcfg.Fs;
t  = fcfg.time;

for iD = 1:numel(fcfg.dataset)
    
    % Settings
    pre   = round(fcfg.pre*fcfg.Fs);
    pos   = round(fcfg.pos*fcfg.Fs);
    
    if size(fcfg.dataset{iD},1) < size(fcfg.dataset{iD},2); fcfg.dataset{iD} = fcfg.dataset{iD}'; end
    var_hld = var(fcfg.dataset{iD});
    [~,thr_lbl] = max(var_hld);
    if size(fcfg.dataset{iD},1) > size(fcfg.dataset{iD},2); fcfg.dataset{iD} = fcfg.dataset{iD}'; end
    
    thresh = max(fcfg.dataset{iD}(thr_lbl,:))/6;
    nsamp = numel(t);
    ntrigchan = 1;
    
    %% Find the timings
    [crsind{1:ntrigchan}] = deal([]);
    for k = 1:ntrigchan
        
        tmpdat  = fcfg.dataset{iD}(k,:);
        ind     = crossing(tmpdat,[],thresh);
        % only keep left edge of pulse
        sel     = (tmpdat(ind+1)-tmpdat(ind-1)) > 0;
        ind     = ind(sel);
        
        crsind{k} = ind;
        clear ind sel tmpdat
    end
    
    bincode = zeros(ntrigchan,nsamp);
    for k = 1:length(crsind); bincode(k,crsind{k}) = 1; end
    detind  = any(bincode,1); % whether there's a detection in any trigchan
    delta   = round(0.01*Fs/2); % make sure detections for the same pulse are grouped before conversion
    ind     = find(detind); % get indices at which a detection occurs in any trigchan
    
    % find detections with a second detection > delta steps before it
    tmp     = cellfun(@(x)any((ind>x-delta & (ind<x))),num2cell(ind));
    R_ind   = find(tmp);
    
    % find detections with a second detection < delta steps after it
    tmp     = cellfun(@(x)any((ind<x+delta & (ind>x))),num2cell(ind));
    L_ind   = find(tmp);
    
    % remove detections that are between two other detections in the same pulse
    Rsel    = ~ismember(R_ind,L_ind);
    Lsel    = ~ismember(L_ind,R_ind);
    R_ind   = R_ind(Rsel); clear Rsel
    L_ind   = L_ind(Lsel); clear Lsel
    
    %% for each pair (L,R), set [ind(L):ind(R)] to 1 in each bincode row w/ a detection
    for k = 1:ntrigchan
        sel = cellfun(@(x,y)any(bincode(k,ind(x):ind(y))),num2cell(L_ind),num2cell(R_ind));
        sel = cellfun(@(x,y)ind(x):ind(y),num2cell(L_ind(sel)),num2cell(R_ind(sel)),'uniformoutput',false);
        sel = unique([sel{:}]);
        bincode(k,sel) = 1;
    end
    
    detind  = any(bincode,1); % whether there's a detection in any trigchan
    % only keep the first pulse in contiguous detections
    samp    = find(detind(1:end-1)==0 & detind(2:end)==1) + 1;
    
    %% Find code for event
    bincode = flipud(bincode(:,samp));
    bincode = mat2cell(bincode,ntrigchan,ones(1,length(samp)));
    % convert binary to decimal
    evt    = cellfun(@(x)polyval(x,2),bincode);
    
    %% Focus on Relevant Events - EJK need to check if works
    if isfield(fcfg,'evt')
        
        ind_evt = zeros(1,numel(evt));
        ind     = [];
        for ievt = 1:length(fcfg.evt)
            ind = [ind find(evt==fcfg.evt(ievt))]; %#ok<AGROW>
        end
        ind_evt(sort(unique(ind))) = 1;
        
        evt(find(~ind_evt)) = []; %#ok<FNDSB>
        samp(find(~ind_evt)) = [];%#ok<FNDSB>
        
    end
    
    %% Remove events too close to the time window
    if isfield(fcfg,'minduration') && ~isempty(samp)
        minsamp = round(fcfg.minduration * Fs);
        sel = find([1 diff(samp)>minsamp]);
        samp = samp(sel);
        evt  = evt(sel);
    elseif isfield(fcfg,'minduration') && ~isempty(samp) && ~isempty(samp)
        samplim = round(fcfg.timelim * Fs);
        sel = find([1 (diff(samp)>samplim(1) & diff(samp)<samplim(2))]);
        samp = samp(sel);
        evt  = evt(sel);
    end
    
    %% Adjust time to make it in the correct sampling frequency
    if isfield(fcfg,'dsfact')
        samp = round(samp/fcfg.dsfact);
        pre = round(pre/fcfg.dsfact);
        pos = round(pos/fcfg.dsfact);
    end
    
    if isfield(fcfg,'ignore') && ~isempty(fcfg.ignore{iD})
        samp(fcfg.ignore{iD}) = [];
        evt(fcfg.ignore{iD})  = [];
    end
    
    %% Add in Events
    if any(strfind(fcfg.trialfun_add{iD},'---'))
        
        trialfun_add_hld = regexp(fcfg.trialfun_add{iD},'---','split');
        
        for iTF = 1:numel(trialfun_add_hld)
            
            evt_hld{iTF} = mmil_readtext([fcfg.prj_dat_hld '/' 'clerical' '/' 'task' '/' fcfg.sbj_nme '/' trialfun_add_hld{iTF}],['\t']);
            if size(evt_hld{iTF},2)==1; evt_hld{iTF} = mmil_readtext([fcfg.prj_dat_hld '/' 'clerical' '/' 'task' '/' fcfg.sbj_nme '/' trialfun_add_hld{iTF}]); end
            evt_hld{iTF}(1,:) = [];
            
        end
        
        evt_hld = cat(1,evt_hld{:});
        
    else
        
        evt_hld = mmil_readtext([fcfg.prj_dat_hld '/' 'clerical' '/' 'task' '/' fcfg.sbj_nme '/' fcfg.trialfun_add{iD}],['\t']);
        if size(evt_hld,2)==1; evt_hld = mmil_readtext([fcfg.prj_dat_hld '/' 'clerical' '/' 'task' '/' fcfg.sbj_nme '/' fcfg.trialfun_add{iD}]); end
        evt_hld(1,:) = [];
        
    end
    
    %    
    evt = cell2mat(evt_hld(:,3));
    
    % Check data
    cell2csv([fcfg.prj_dat_hld '/' 'clerical' '/' 'task' '/' fcfg.sbj_nme '/' 'TimeTest' '_' fcfg.trialfun_add{iD}],num2cell([diff(samp)' diff(cell2mat(evt_hld(:,2))) diff(samp)'-diff(cell2mat(evt_hld(:,2)))  ]));
         
    %% Make new trl file
    trl{iD} = zeros(length(evt),4);
    for i=1:length(evt)
        trl{iD}(i,1) = samp(i)+pre;
        trl{iD}(i,2) = samp(i)+pos;
        trl{iD}(i,3) = pre;
        trl{iD}(i,4) = evt(i);
    end
end

end
%% Test ability to find correct triggers
% figure()
% hold on
% plot(dat.trial{1}(1,:)/3463134.48017823 + 14)
% ylim([0 16]);
% plot(dat.trial{1}(2,:)/3463134.48017823 + 12)
% plot(dat.trial{1}(3,:)/3463134.48017823 + 10)
% plot(dat.trial{1}(4,:)/3463134.48017823 + 8)
% plot(dat.trial{1}(5,:)/3463134.48017823 + 6)
% plot(dat.trial{1}(6,:)/3463134.48017823 + 4)
% plot(dat.trial{1}(7,:)/3463134.48017823 + 2)
% plot(dat.trial{1}(8,:)/3463134.48017823 + 0)
% 
% for it = 1:length(samp)
%     vline(a(it),'g');
% end


function FW_events(fcfg)

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.sbj_dat_hld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

end_dir     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'end_dir'); %mmil_readtext([fcfg.clr_fld '/end_dir/' subj]);
tsk         = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'tsk'); %mmil_readtext([fcfg.clr_fld '/tsk/' subj]);     

if strcmpi(end_dir,'eeg') | strcmpi(end_dir,'edf')
    eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes.csv']);
elseif strcmpi(end_dir,'mat') | strcmpi(end_dir,'.mat')
    eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes.csv']);
end

eve_fwv = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events.csv']);
if strcmpi(tsk,'fw_1')
    idn_col = 2;
    eve_col = 4;
else
    idn_col = 4;
    eve_col = 3;
end

%% Check data
wrg_num = 1;
wrg_mtc = 1;

while wrg_num || wrg_mtc
    
    if ~exist([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv'])
        if strcmpi(end_dir,'eeg') | strcmpi(end_dir,'edf')
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes_corrected.csv'],{});
        elseif strcmpi(end_dir,'mat') | strcmpi(end_dir,'.mat')
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes_corrected.csv'],{});
        end
        cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv'],{})
    else
        if strcmpi(end_dir,'eeg') | strcmpi(end_dir,'edf')
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes_corrected.csv']);
        elseif strcmpi(end_dir,'mat') | strcmpi(end_dir,'.mat')
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes_corrected.csv']);
        end
        eve_fwv = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv']);
    end
        
    if size(eve_cde,1) ~= size(eve_fwv,1)
        disp([fcfg.sbj_nme ': number of trials do no match' ])
        pause()
        pause(5)       
    elseif ~all(cell2mat(eve_cde(:,1)) == cell2mat(eve_fwv(:,4)))
        wrg_num = 0;
        
        disp([fcfg.sbj_nme ': events do no match' ])
        find(cell2mat(eve_cde(:,1)) ~= cell2mat(eve_fwv(:,4)))
        pause()
        pause(5)
    else
        wrg_num = 0;
        wrg_mtc = 0;
    end

    
    try isempty(mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv'])); catch; emp_txt = 1; end
    if ~wrg_num && ~wrg_mtc && isvar('emp_txt') && emp_txt
        if strcmpi(end_dir,'eeg') | strcmpi(end_dir,'edf')
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes_corrected.csv'],eve_cde);
        elseif strcmpi(end_dir,'mat') | strcmpi(end_dir,'.mat')
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes_corrected.csv'],eve_cde);
        end
        cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv'],eve_fwv)
    else
        if strcmpi(end_dir,'eeg') | strcmpi(end_dir,'edf')
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes_corrected.csv']);
        elseif strcmpi(end_dir,'mat') | strcmpi(end_dir,'.mat')
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes_corrected.csv']);
        end
        eve_fwv = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv']);
    end
    
end

%% Remove missing trials & Reorder trials
for iD = 1:numel(bcc_dat.data_name)
    if any(bcc_dat.(bcc_dat.data_name{iD}).trialinfo==4) && any(bcc_dat.(bcc_dat.data_name{iD}).trialinfo==8) && any(bcc_dat.(bcc_dat.data_name{iD}).trialinfo==16) && any(bcc_dat.(bcc_dat.data_name{iD}).trialinfo==32)
        for iFW = 1:size(bcc_dat.(bcc_dat.data_name{iD}).trialinfo,1)
            bcc_dat.(bcc_dat.data_name{iD}).trialinfo(iFW) = numel(factor(bcc_dat.(bcc_dat.data_name{iD}).trialinfo(iFW)))+1;
        end
    end
end
    
if strcmpi(end_dir,'eeg') | strcmpi(end_dir,'edf')
    eve_cde_org = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes.csv']);
elseif strcmpi(end_dir,'mat') | strcmpi(end_dir,'.mat')
    eve_cde_org = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes.csv']);
end

if sum(bcc_dat.(bcc_dat.data_name{1}).trialinfo==3) > sum(bcc_dat.(bcc_dat.data_name{1}).trialinfo==8)
    inc_trl = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==3 | bcc_dat.(bcc_dat.data_name{1}).trialinfo==4 | bcc_dat.(bcc_dat.data_name{1}).trialinfo==5 | bcc_dat.(bcc_dat.data_name{1}).trialinfo==6);
elseif sum(bcc_dat.(bcc_dat.data_name{1}).trialinfo==3) < sum(bcc_dat.(bcc_dat.data_name{1}).trialinfo==8)
    inc_trl = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==8 | bcc_dat.(bcc_dat.data_name{1}).trialinfo==4 | bcc_dat.(bcc_dat.data_name{1}).trialinfo==5 | bcc_dat.(bcc_dat.data_name{1}).trialinfo==6);
end

cfg = [];
cfg.trials = inc_trl;
bcc_dat = ft_func(@ft_preprocessing,cfg,bcc_dat);

if sum(bcc_dat.(bcc_dat.data_name{1}).trialinfo==3) < sum(bcc_dat.(bcc_dat.data_name{1}).trialinfo==8)
    bcc_dat.(bcc_dat.data_name{1}).trialinfo(bcc_dat.(bcc_dat.data_name{1}).trialinfo==8) = 3;
end
   
thr_lst_ind = cell2mat(eve_fwv(cell2mat(eve_fwv(:,4))==3,3)); thr_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==3,2)); thr_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==3,2)); thr_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==3);
for_lst_ind = cell2mat(eve_fwv(cell2mat(eve_fwv(:,4))==4,3)); for_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==4,2)); for_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==4,2)); for_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==4);
fve_lst_ind = cell2mat(eve_fwv(cell2mat(eve_fwv(:,4))==5,3)); fve_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==5,2)); fve_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==5,2)); fve_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==5);
six_lst_ind = cell2mat(eve_fwv(cell2mat(eve_fwv(:,4))==6,3)); six_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==6,2)); six_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==6,2)); six_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==6);

[~,thr_ext] = intersect(thr_cde_org_ind,setdiff(thr_cde_org_ind,thr_cde_ind)); if ~isempty(thr_ext) && thr_ext(end) <= numel(thr_trl_ind); thr_ext = thr_trl_ind(thr_ext); else thr_ext = []; end
[~,for_ext] = intersect(for_cde_org_ind,setdiff(for_cde_org_ind,for_cde_ind)); if ~isempty(for_ext) && for_ext(end) <= numel(for_trl_ind); for_ext = for_trl_ind(for_ext); else for_ext = []; end
[~,fve_ext] = intersect(fve_cde_org_ind,setdiff(fve_cde_org_ind,fve_cde_ind)); if ~isempty(fve_ext) && fve_ext(end) <= numel(fve_trl_ind); fve_ext = fve_trl_ind(fve_ext); else fve_ext = []; end
[~,six_ext] = intersect(six_cde_org_ind,setdiff(six_cde_org_ind,six_cde_ind)); if ~isempty(six_ext) && six_ext(end) <= numel(six_trl_ind); six_ext = six_trl_ind(six_ext); else six_ext = []; end

cfg = [];
cfg.trials = 1:numel(bcc_dat.(bcc_dat.data_name{1}).trialinfo);
cfg.trials([thr_ext ; for_ext ; fve_ext ; six_ext]) = [];
bcc_dat = ft_func(@ft_preprocessing,cfg,bcc_dat);

thr_ind = 1; for_ind = 1; fve_ind = 1; six_ind = 1; eve_key = zeros(numel(bcc_dat.(bcc_dat.data_name{1}).trialinfo),1);

rmv_ind = [];
for iT = 1:numel(bcc_dat.(bcc_dat.data_name{1}).trialinfo)
    
    switch bcc_dat.(bcc_dat.data_name{1}).trialinfo(iT)
        case 3
            try
                eve_key(iT) = thr_lst_ind(thr_ind);
                thr_ind = thr_ind + 1;
            catch
                rmv_ind = [rmv_ind iT];
            end
        case 4
            try
                eve_key(iT) = for_lst_ind(for_ind);
                for_ind = for_ind + 1;
            catch
                rmv_ind = [rmv_ind iT];
            end
        case 5
            try
                eve_key(iT) = fve_lst_ind(fve_ind);
                fve_ind = fve_ind + 1;
            catch
                rmv_ind = [rmv_ind iT];
            end
        case 6
            try
                eve_key(iT) = six_lst_ind(six_ind);
                six_ind = six_ind + 1;
            catch
                rmv_ind = [rmv_ind iT];
            end
    end
end

%% Add in events
stt = mmil_readtext([fcfg.clr_fld '/' 'task' '/' tsk{1} '.csv']);

for iD = 1:numel(bcc_dat.data_name)
    
    bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo = bcc_dat.(bcc_dat.data_name{1}).trialinfo;
    
    bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.vis_ovr = ones(numel(bcc_dat.(bcc_dat.data_name{1}).trialinfo),1)*101;
    
    bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.stm_idn = stt(eve_key,1);
    
    bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.wrd_lng = cell2mat(stt(eve_key,2));
        bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.wrd_lng(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo~=3) = nan;
    bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.wrd_frq = cell2mat(stt(eve_key,3));
        bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.wrd_frq(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo~=3) = nan;
    bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.wrd_ngh = cell2mat(stt(eve_key,4));
        bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.wrd_ngh(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo~=3) = nan;
    bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.wrd_big = cell2mat(stt(eve_key,5));
        bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.wrd_big(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo~=3) = nan;
        
    bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.nwn_lng = cell2mat(stt(eve_key,2));
        bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.nwn_lng(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo~=5) = nan;
    bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.nwn_ngh = cell2mat(stt(eve_key,4));
        bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.nwn_ngh(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo~=5) = nan;
    bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.nwn_big = cell2mat(stt(eve_key,5));
        bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.nwn_big(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo~=5) = nan;
    
end

%% Make additional events
for iD = 2:numel(bcc_dat.data_name); bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve; end

cfg = [];
cfg.eve_rul     = {'med_spl'           'med_spl'           'med_spl'           'med_spl'};
cfg.new_events  = {[111 111 0 112 112] [121 121 0 122 122] [131 132]           [141 141 0 142 142]};
cfg.crt_alt_eve = {'wrd_lng_med'       'wrd_frq_med'       'wrd_ngh_med'       'wrd_big_med'};
cfg.fld_nme     = {'wrd_lng'           'wrd_frq'           'wrd_ngh'           'wrd_big'};
bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);
    
cfg = [];
cfg.eve_rul     = {'med_spl'           'med_spl'           };
cfg.new_events  = {[211 211 0 212 212] [221 221 0 222 222] };
cfg.crt_alt_eve = {'nwn_lng_med'       'nwn_big_med'       };
cfg.fld_nme     = {'nwn_lng'           'nwn_big'           };
bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);    

%% Save Stats
if ~exist([fcfg.clr_fld '/' 'events'],'dir'); mkdir([fcfg.clr_fld '/' 'events']); end
if ~exist([fcfg.clr_fld '/' 'events' '/' fcfg.sbj_nme],'dir'); mkdir([fcfg.clr_fld '/' 'events' '/' fcfg.sbj_nme]); end

eve_hld = fieldnames(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve);
eve_hld = sort(eve_hld);
eve_tbl = num2cell(eve_key);
for iEH = 1:numel(eve_hld)
    if isnumeric(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_hld{iEH})(1))
        eve_tbl(:,iEH+1) = num2cell(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_hld{iEH}));
    else
        eve_tbl(:,iEH+1) = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_hld{iEH});
    end
end

cell2csv([fcfg.clr_fld '/' 'events' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '.csv'],[[{''} eve_hld'] ; eve_tbl])

%% Save Data
cfg = [];
cfg.str_nme  = 'fwv_dat';
cfg.save     = 'yes';
cfg.filename =[fcfg.sbj_dat_hld '/' fcfg.sbj_nme '_overall_data.mat'];
ft_func([],cfg,bcc_dat);

end






function iSASZ_initial_events3(fcfg)

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.sbj_dat_hld '/' 'epoch_data' '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

end_dir     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'end_dir'); %mmil_readtext([fcfg.clr_fld '/end_dir/' subj]);

% Load codes
if ~isempty(string_find(end_dir,'eeg')) | ~isempty(string_find(end_dir,'set'))
    eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes.csv']);
elseif ~isempty(string_find(end_dir,'mat'))
    eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes.csv']);
elseif ~isempty(string_find(end_dir,'edf'))
    eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_edf_codes.csv']);
end

% Load real events
stm_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events.csv']);
eve_col = 3;

%% Check data
wrg_num = 1;
wrg_mtc = 1;

while wrg_num || wrg_mtc
    
    if ~exist([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv'])
        if  ~isempty(string_find(end_dir,'eeg')) | ~isempty(string_find(end_dir,'set'))
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes_corrected.csv'],{});
        elseif ~isempty(string_find(end_dir,'mat'))
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes_corrected.csv'],{});
        elseif ~isempty(string_find(end_dir,'edf'))
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_edf_codes_corrected.csv'],{});
        end
        cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv'],{})
    else
        if  ~isempty(string_find(end_dir,'eeg')) | ~isempty(string_find(end_dir,'set'))
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes_corrected.csv']);
        elseif ~isempty(string_find(end_dir,'mat'))
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes_corrected.csv']);
        elseif ~isempty(string_find(end_dir,'edf'))
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_edf_codes_corrected.csv']);
        end
        stm_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv']);
    end
    
    if size(eve_cde,1) ~= size(stm_cde,1)
        disp([fcfg.sbj_nme ': number of trials do no match' ])
        pause()
        pause(5)
    elseif ~all(cell2mat(eve_cde(:,1)) == cell2mat(stm_cde(:,eve_col)))
        wrg_num = 0;
        
        disp([fcfg.sbj_nme ': events do no match' ])
        find(cell2mat(eve_cde(:,1)) ~= cell2mat(stm_cde(:,eve_col)))
        pause()
        pause(5)
    else
        wrg_num = 0;
        wrg_mtc = 0;
    end
    
    
    try isempty(mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv'])); catch; emp_txt = 1; end
    if ~wrg_num && ~wrg_mtc && isvar('emp_txt') && emp_txt
        if ~isempty(string_find(end_dir,'eeg')) | ~isempty(string_find(end_dir,'set'))
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes_corrected.csv'],eve_cde);
        elseif ~isempty(string_find(end_dir,'mat'))
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes_corrected.csv'],eve_cde);
            elseif ~isempty(string_find(end_dir,'edf'))
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_edf_codes_corrected.csv'],eve_cde);
        end
        cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv'],stm_cde)
    else
        if ~isempty(string_find(end_dir,'eeg')) | ~isempty(string_find(end_dir,'set'))
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes_corrected.csv']);
        elseif ~isempty(string_find(end_dir,'mat'))
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes_corrected.csv']);
            elseif ~isempty(string_find(end_dir,'edf'))
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_edf_codes_corrected.csv']);
        end
        stm_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv']);
    end
    
end

%% Remove missing trials & Reorder trials
if ~isempty(string_find(end_dir,'eeg'))  | ~isempty(string_find(end_dir,'set'))
    eve_cde_org = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes.csv']);
elseif ~isempty(string_find(end_dir,'mat'))
    eve_cde_org = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes.csv']);
elseif ~isempty(string_find(end_dir,'edf'))
    eve_cde_org = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_edf_codes.csv']);
end

eve_key     = [ 11  10  13  12  35   95   115  54    104  27   47   76   86   66     ...
    112 111 122 121 8212 8222 5212 5222 3211 3221  10211 10221 9211 9221 1212 1222 7222 2222 4221 6221 6211; ...
    1   2   3   4   5    5    5    6    6    7     7     8     8    8          ...
    11  12  13  14  15   15   15   15   16   16    16    16    16   16   17   17 17   17   18   18   18    ];

for iL = 1:numel(bcc_dat.data_name)
    if any(bcc_dat.(bcc_dat.data_name{iL}).trialinfo>20)
        if any(bcc_dat.(bcc_dat.data_name{iL}).trialinfo==54) || any(bcc_dat.(bcc_dat.data_name{iL}).trialinfo==121)
            for iE = 1:size(eve_key,2)
                bcc_dat.(bcc_dat.data_name{iL}).trialinfo(bcc_dat.(bcc_dat.data_name{iL}).trialinfo==eve_key(1,iE)) = eve_key(2,iE);
            end
        end
    end
end

szv_one_lst_ind = cell2mat(stm_cde(cell2mat(stm_cde(:,3))==1,1));  szv_one_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==1,2));  szv_one_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==1,2));  szv_one_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==1);
szv_two_lst_ind = cell2mat(stm_cde(cell2mat(stm_cde(:,3))==2,1));  szv_two_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==2,2));  szv_two_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==2,2));  szv_two_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==2);
szv_thr_lst_ind = cell2mat(stm_cde(cell2mat(stm_cde(:,3))==3,1));  szv_thr_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==3,2));  szv_thr_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==3,2));  szv_thr_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==3);
szv_for_lst_ind = cell2mat(stm_cde(cell2mat(stm_cde(:,3))==4,1));  szv_for_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==4,2));  szv_for_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==4,2));  szv_for_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==4);
szv_fve_lst_ind = cell2mat(stm_cde(cell2mat(stm_cde(:,3))==5,1));  szv_fve_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==5,2));  szv_fve_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==5,2));  szv_fve_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==5);
szv_six_lst_ind = cell2mat(stm_cde(cell2mat(stm_cde(:,3))==6,1));  szv_six_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==6,2));  szv_six_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==6,2));  szv_six_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==6);
szv_svn_lst_ind = cell2mat(stm_cde(cell2mat(stm_cde(:,3))==7,1));  szv_svn_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==7,2));  szv_svn_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==7,2));  szv_svn_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==7);
szv_eig_lst_ind = cell2mat(stm_cde(cell2mat(stm_cde(:,3))==8,1));  szv_eig_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==8,2));  szv_eig_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==8,2));  szv_eig_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==8);

sza_one_lst_ind = cell2mat(stm_cde(cell2mat(stm_cde(:,3))==11,1)); sza_one_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==11,2)); sza_one_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==11,2)); sza_one_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==11);
sza_two_lst_ind = cell2mat(stm_cde(cell2mat(stm_cde(:,3))==12,1)); sza_two_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==12,2)); sza_two_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==12,2)); sza_two_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==12);
sza_thr_lst_ind = cell2mat(stm_cde(cell2mat(stm_cde(:,3))==13,1)); sza_thr_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==13,2)); sza_thr_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==13,2)); sza_thr_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==13);
sza_for_lst_ind = cell2mat(stm_cde(cell2mat(stm_cde(:,3))==14,1)); sza_for_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==14,2)); sza_for_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==14,2)); sza_for_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==14);
sza_fve_lst_ind = cell2mat(stm_cde(cell2mat(stm_cde(:,3))==15,1)); sza_fve_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==15,2)); sza_fve_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==15,2)); sza_fve_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==15);
sza_six_lst_ind = cell2mat(stm_cde(cell2mat(stm_cde(:,3))==16,1)); sza_six_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==16,2)); sza_six_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==16,2)); sza_six_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==16);
sza_svn_lst_ind = cell2mat(stm_cde(cell2mat(stm_cde(:,3))==17,1)); sza_svn_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==17,2)); sza_svn_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==17,2)); sza_svn_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==17);
sza_eig_lst_ind = cell2mat(stm_cde(cell2mat(stm_cde(:,3))==18,1)); sza_eig_cde_ind = cell2mat(eve_cde(cell2mat(eve_cde(:,1))==18,2)); sza_eig_cde_org_ind = cell2mat(eve_cde_org(cell2mat(eve_cde_org(:,1))==18,2)); sza_eig_trl_ind = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==18);

[~,szv_one_ext] = intersect(szv_one_cde_org_ind,setdiff(szv_one_cde_org_ind,szv_one_cde_ind)); if ~isempty(szv_one_ext) && szv_one_ext(end) <= numel(szv_one_trl_ind); szv_one_ext = szv_one_trl_ind(szv_one_ext); else szv_one_ext = []; end
[~,szv_two_ext] = intersect(szv_two_cde_org_ind,setdiff(szv_two_cde_org_ind,szv_two_cde_ind)); if ~isempty(szv_two_ext) && szv_two_ext(end) <= numel(szv_two_trl_ind); szv_two_ext = szv_two_trl_ind(szv_two_ext); else szv_two_ext = []; end
[~,szv_thr_ext] = intersect(szv_thr_cde_org_ind,setdiff(szv_thr_cde_org_ind,szv_thr_cde_ind)); if ~isempty(szv_thr_ext) && szv_thr_ext(end) <= numel(szv_thr_trl_ind); szv_thr_ext = szv_thr_trl_ind(szv_thr_ext); else szv_thr_ext = []; end
[~,szv_for_ext] = intersect(szv_for_cde_org_ind,setdiff(szv_for_cde_org_ind,szv_for_cde_ind)); if ~isempty(szv_for_ext) && szv_for_ext(end) <= numel(szv_for_trl_ind); szv_for_ext = szv_for_trl_ind(szv_for_ext); else szv_for_ext = []; end
[~,szv_fve_ext] = intersect(szv_fve_cde_org_ind,setdiff(szv_fve_cde_org_ind,szv_fve_cde_ind)); if ~isempty(szv_fve_ext) && szv_fve_ext(end) <= numel(szv_fve_trl_ind); szv_fve_ext = szv_fve_trl_ind(szv_fve_ext); else szv_fve_ext = []; end
[~,szv_six_ext] = intersect(szv_six_cde_org_ind,setdiff(szv_six_cde_org_ind,szv_six_cde_ind)); if ~isempty(szv_six_ext) && szv_six_ext(end) <= numel(szv_six_trl_ind); szv_six_ext = szv_six_trl_ind(szv_six_ext); else szv_six_ext = []; end
[~,szv_svn_ext] = intersect(szv_svn_cde_org_ind,setdiff(szv_svn_cde_org_ind,szv_svn_cde_ind)); if ~isempty(szv_svn_ext) && szv_svn_ext(end) <= numel(szv_svn_trl_ind); szv_svn_ext = szv_svn_trl_ind(szv_svn_ext); else szv_svn_ext = []; end
[~,szv_eig_ext] = intersect(szv_eig_cde_org_ind,setdiff(szv_eig_cde_org_ind,szv_eig_cde_ind)); if ~isempty(szv_eig_ext) && szv_eig_ext(end) <= numel(szv_eig_trl_ind); szv_eig_ext = szv_eig_trl_ind(szv_eig_ext); else szv_eig_ext = []; end

[~,sza_one_ext] = intersect(sza_one_cde_org_ind,setdiff(sza_one_cde_org_ind,sza_one_cde_ind)); if ~isempty(sza_one_ext) && sza_one_ext(end) <= numel(sza_one_trl_ind); sza_one_ext = sza_one_trl_ind(sza_one_ext); else sza_one_ext = []; end
[~,sza_two_ext] = intersect(sza_two_cde_org_ind,setdiff(sza_two_cde_org_ind,sza_two_cde_ind)); if ~isempty(sza_two_ext) && sza_two_ext(end) <= numel(sza_two_trl_ind); sza_two_ext = sza_two_trl_ind(sza_two_ext); else sza_two_ext = []; end
[~,sza_thr_ext] = intersect(sza_thr_cde_org_ind,setdiff(sza_thr_cde_org_ind,sza_thr_cde_ind)); if ~isempty(sza_thr_ext) && sza_thr_ext(end) <= numel(sza_thr_trl_ind); sza_thr_ext = sza_thr_trl_ind(sza_thr_ext); else sza_thr_ext = []; end
[~,sza_for_ext] = intersect(sza_for_cde_org_ind,setdiff(sza_for_cde_org_ind,sza_for_cde_ind)); if ~isempty(sza_for_ext) && sza_for_ext(end) <= numel(sza_for_trl_ind); sza_for_ext = sza_for_trl_ind(sza_for_ext); else sza_for_ext = []; end
[~,sza_fve_ext] = intersect(sza_fve_cde_org_ind,setdiff(sza_fve_cde_org_ind,sza_fve_cde_ind)); if ~isempty(sza_fve_ext) && sza_fve_ext(end) <= numel(sza_fve_trl_ind); sza_fve_ext = sza_fve_trl_ind(sza_fve_ext); else sza_fve_ext = []; end
[~,sza_six_ext] = intersect(sza_six_cde_org_ind,setdiff(sza_six_cde_org_ind,sza_six_cde_ind)); if ~isempty(sza_six_ext) && sza_six_ext(end) <= numel(sza_six_trl_ind); sza_six_ext = sza_six_trl_ind(sza_six_ext); else sza_six_ext = []; end
[~,sza_svn_ext] = intersect(sza_svn_cde_org_ind,setdiff(sza_svn_cde_org_ind,sza_svn_cde_ind)); if ~isempty(sza_svn_ext) && sza_svn_ext(end) <= numel(sza_svn_trl_ind); sza_svn_ext = sza_svn_trl_ind(sza_svn_ext); else sza_svn_ext = []; end
[~,sza_eig_ext] = intersect(sza_eig_cde_org_ind,setdiff(sza_eig_cde_org_ind,sza_eig_cde_ind)); if ~isempty(sza_eig_ext) && sza_eig_ext(end) <= numel(sza_eig_trl_ind); sza_eig_ext = sza_eig_trl_ind(sza_eig_ext); else sza_eig_ext = []; end

cfg = [];
cfg.trials = 1:numel(bcc_dat.(bcc_dat.data_name{1}).trialinfo);
cfg.trials([ szv_one_ext ; szv_two_ext ; szv_thr_ext ; szv_for_ext ; szv_fve_ext ; szv_six_ext ; szv_svn_ext ; szv_eig_ext ; ...
    sza_one_ext ; sza_two_ext ; sza_thr_ext ; sza_for_ext ; sza_fve_ext ; sza_six_ext ; sza_svn_ext ; sza_eig_ext ]) = [];
bcc_dat = ft_func(@ft_preprocessing,cfg,bcc_dat);

szv_one_ind = 1; szv_two_ind = 1; szv_thr_ind = 1; szv_for_ind = 1; szv_fve_ind = 1; szv_six_ind = 1; szv_svn_ind = 1; szv_eig_ind = 1;
sza_one_ind = 1; sza_two_ind = 1; sza_thr_ind = 1; sza_for_ind = 1; sza_fve_ind = 1; sza_six_ind = 1; sza_svn_ind = 1; sza_eig_ind = 1;
eve_key = zeros(numel(bcc_dat.(bcc_dat.data_name{1}).trialinfo),1);

for iT = 1:numel(bcc_dat.(bcc_dat.data_name{1}).trialinfo)
    
    switch bcc_dat.(bcc_dat.data_name{1}).trialinfo(iT)
        case 1
            eve_key(iT) = szv_one_lst_ind(szv_one_ind);
            szv_one_ind = szv_one_ind + 1;
        case 2
            eve_key(iT) = szv_two_lst_ind(szv_two_ind);
            szv_two_ind = szv_two_ind + 1;
        case 3
            eve_key(iT) = szv_thr_lst_ind(szv_thr_ind);
            szv_thr_ind = szv_thr_ind + 1;
        case 4
            eve_key(iT) = szv_for_lst_ind(szv_for_ind);
            szv_for_ind = szv_for_ind + 1;
        case 5
            eve_key(iT) = szv_fve_lst_ind(szv_fve_ind);
            szv_fve_ind = szv_fve_ind + 1;
        case 6
            eve_key(iT) = szv_six_lst_ind(szv_six_ind);
            szv_six_ind = szv_six_ind + 1;
        case 7
            eve_key(iT) = szv_svn_lst_ind(szv_svn_ind);
            szv_svn_ind = szv_svn_ind + 1;
        case 8
            eve_key(iT) = szv_eig_lst_ind(szv_eig_ind);
            szv_eig_ind = szv_eig_ind + 1;
        case 11
            eve_key(iT) = sza_one_lst_ind(sza_one_ind);
            sza_one_ind = sza_one_ind + 1;
        case 12
            eve_key(iT) = sza_two_lst_ind(sza_two_ind);
            sza_two_ind = sza_two_ind + 1;
        case 13
            eve_key(iT) = sza_thr_lst_ind(sza_thr_ind);
            sza_thr_ind = sza_thr_ind + 1;
        case 14
            eve_key(iT) = sza_for_lst_ind(sza_for_ind);
            sza_for_ind = sza_for_ind + 1;
        case 15
            eve_key(iT) = sza_fve_lst_ind(sza_fve_ind);
            sza_fve_ind = sza_fve_ind + 1;
        case 16
            eve_key(iT) = sza_six_lst_ind(sza_six_ind);
            sza_six_ind = sza_six_ind + 1;
        case 17
            eve_key(iT) = sza_svn_lst_ind(sza_svn_ind);
            sza_svn_ind = sza_svn_ind + 1;
        case 18
            eve_key(iT) = sza_eig_lst_ind(sza_eig_ind);
            sza_eig_ind = sza_eig_ind + 1;
    end
end

if any(eve_key==0)
    eve_key
    disp('eve_key wrong, check it out')
    pause()
    eve_key = cell2mat(stm_cde(:,1));
end

%% Add in events
stt = mmil_readtext([fcfg.clr_fld '/' 'task' '/' 'tsk.csv']);

for iD = 1:numel(bcc_dat.data_name)
    
    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo = cell2mat(stt(eve_key,3));
    
    % Add in identity
    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.stm_idn = stt(eve_key,2);
    
    % Change initial repeated words
    rep_wrd = tabulate(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.stm_idn);
    rep_wrd = rep_wrd(cell2mat(rep_wrd(:,2))>5,1);
    
    for iR = 1:numel(rep_wrd)
        
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo(find(strcmp(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.stm_idn,rep_wrd{iR}),1)) = ...
            bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo(find(strcmp(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.stm_idn,rep_wrd{iR}),1))-4;
    end
    
    % Add in orthographic statistics
    if any(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo < 10)
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_lng = cell2mat(stt(eve_key,4)); % vis length
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_lng(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo>10) = NaN;
        
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_big = cell2mat(stt(eve_key,7)); % vis bigram
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_big(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo>10) = NaN;
        
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_frq = cell2mat(stt(eve_key,5)); % vis freq
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_frq(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo>10) = NaN;
        
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_ngh = cell2mat(stt(eve_key,6)); % vis ngh
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_ngh(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo>10) = NaN;
        
        % Add in semantic priming
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_sem_prm = cell2mat(stt(eve_key,8));
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_sem_prm(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo>10) = NaN;
    end
    
    % Add in phonetic statistics
    if any(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo > 10)
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_lng = cell2mat(stt(eve_key,4)); % aud length
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_lng(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo<10) = NaN;
        
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_big = cell2mat(stt(eve_key,7)); % aud bigram
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_big(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo<10) = NaN;
        
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_frq = cell2mat(stt(eve_key,5)); % aud freq
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_frq(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo<10) = NaN;
        
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_ngh = cell2mat(stt(eve_key,6)); % aud ngh
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_ngh(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo<10) = NaN;
        
        % Add in semantic priming       
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_sem_prm = cell2mat(stt(eve_key,8));
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_sem_prm(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo<10) = NaN;
    end
    
end

%% Update Events
if any(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo < 10)
    
    % Overall
    for iD = 1:numel(bcc_dat.data_name)
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_ovr = (bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo<10) + 100;
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_ovr(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_ovr==100) = 0;
    end
    
    % Old/New & Ani/Obj
    cfg = [];
    cfg.eve_rul     = { 'new'                 'new'         };
    cfg.old_events  = { {[1 2 3 4] [5 6 7 8]} {[1 3] [2 4]} };
    cfg.new_events  = { [111 112]             [121 122]     };
    cfg.crt_alt_eve = { 'vis_new_old'         'vis_ani_obj' };
    cfg.use_alt_eve = { 'trialinfo'           'trialinfo'   };
    bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);
    
    % Add in semantic priming
    cfg = [];
    cfg.eve_rul     = {'med_spl'           };
    cfg.new_events  = {[123 123 0 124 124] };
    cfg.crt_alt_eve = {'vis_sem_prm_med'   };
    cfg.fld_nme     = {'vis_sem_prm'       };
    bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);
    
    % Add in orthographic statistics
    cfg = [];
    cfg.eve_rul     = {'med_spl'           'med_spl'           'med_spl'           'med_spl'};
    cfg.new_events  = {[131 131 0 132 132] [141 141 0 142 142] [151 151 0 152 152] [161 161 162 162]};
    cfg.crt_alt_eve = {'vis_lng_med'       'vis_big_med'       'vis_frq_med'       'vis_ngh_med'};
    cfg.fld_nme     = {'vis_lng'           'vis_big'           'vis_frq'           'vis_ngh'};
    bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);
    
end

if any(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo > 10)
    
    % Overall
    for iD = 1:numel(bcc_dat.data_name)
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_ovr = (bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo>10) + 200;
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_ovr(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_ovr==200) = 0;
    end
    
    % Old/New & Ani/Obj
    cfg = [];
    cfg.eve_rul     = { 'new'                         'new'};
    cfg.old_events  = { {[11 12 13 14] [15 16 17 18]} {[11 13] [12 14]} };
    cfg.new_events  = { [211 212]                     [221 222]};
    cfg.crt_alt_eve = { 'aud_new_old'                 'aud_ani_obj'};
    cfg.use_alt_eve = { 'trialinfo'                   'trialinfo'};
    bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);
    
    % Add in semantic priming
    cfg = [];
    cfg.eve_rul     = {'med_spl'};
    cfg.new_events  = {[223 223 0 224 224]};
    cfg.crt_alt_eve = {'aud_sem_prm_med' };
    cfg.fld_nme     = {'aud_sem_prm' };
    bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);
    
    % Add in phonetic statistics
    cfg = [];
    cfg.eve_rul     = {'med_spl'           'med_spl'           'med_spl'           'med_spl' };
    cfg.new_events  = {[231 231 0 232 232] [241 241 0 242 242] [251 251 0 252 252] [261 261 0 262 262]};
    cfg.crt_alt_eve = {'aud_lng_med'       'aud_big_med'       'aud_frq_med'       'aud_ngh_med'};
    cfg.fld_nme     = {'aud_lng'           'aud_big'           'aud_frq'           'aud_ngh'};
    bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);
    
end

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
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.filename = [fcfg.sbj_dat_hld '/' 'epoch_data' '/' fcfg.sbj_nme '_overall_data.mat'];
ft_func([],cfg,bcc_dat);

end
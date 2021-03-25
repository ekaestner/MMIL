function sl_events(fcfg)

%% Initial Load
sbj  = fcfg.sbj_nme;

fprintf([fcfg.sbj_nme ': Starting initial stats work on %s \n'],sbj)

infile = [fcfg.out_pth sbj '_overall_data.mat'];
outpath = fcfg.out_pth;

cfg = [];
cfg.load = 'yes';
cfg.file = [outpath '/' sbj '_overall_data.mat'];
bcc_dat  = ft_func([],cfg);

tsk         = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'tsk'); %mmil_readtext([fcfg.clr_fld '/tsk/' subj]);     

try
    
    sbj_lst = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_list.csv']);
    sbj_log = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_log.csv']);

catch
    
    cfg = [];
    cfg.typ = 'file';
    lst_fle = mmil_find_file(cfg,[fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' 'list' '/' '*']);
    sbj_lst = [];
    for iL = 1:numel(lst_fle)
        sbj_lst = [sbj_lst ; mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' 'list' '/' lst_fle{iL}])];
    end
    sbj_lst = [num2cell(1:size(sbj_lst,1))' sbj_lst];
    
    cfg = [];
    cfg.typ = 'file';
    log_fle = mmil_find_file(cfg,[fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' 'log' '/' '*']);
    sbj_log = [];
    for iL = 1:numel(log_fle)
        
        log_hld = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' 'log' '/' log_fle{iL}],'[\t]');
        log_hld = log_hld(6:end,[3 4 ]);
        rmv_ind = [];
        for iT = 1:size(log_hld,1)
            if isempty(log_hld{iT,2}) || isstr(log_hld{iT,2}) || ~(log_hld{iT,2}==1 || log_hld{iT,2}==2 || log_hld{iT,2}==3 || log_hld{iT,2}==4)
                rmv_ind = [rmv_ind ; iT];
            end
        end
        log_hld(rmv_ind,:) = [];
        
        sbj_log = [sbj_log ; log_hld];
    end
    sbj_log = sbj_log(:,2);
    
    cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_list.csv'],sbj_lst);
    cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_log.csv'],sbj_log);
end

if size(sbj_lst,1) ~= size(sbj_log,1)
    disp(['FIX EVENT LENGTH FOR ' fcfg.sbj_nme])
    pause()
    sbj_lst = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_list_cor.csv']);
    sbj_log = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_log_cor.csv']);
end

if ~all(cell2mat(sbj_lst(:,4)) == cell2mat(sbj_log(:,1)))
    disp(['FIX EVENT NUMBERS FOR ' fcfg.sbj_nme])
    pause()
    sbj_lst = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_list_cor.csv']);
    sbj_log = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_log_cor.csv']);
end

sbj_eve_ind = cell2mat(sbj_lst(:,1));

%% Main Event update
for iE = 1:numel(bcc_dat.data_name)
    if size(sbj_eve_ind,1) == size(bcc_dat.(bcc_dat.data_name{iE}).trialinfo,1)
        
        bcc_dat.(bcc_dat.data_name{iE}).trialinfo = cell2mat(sbj_lst(sbj_eve_ind,4));
        bcc_dat.(bcc_dat.data_name{iE}).cfg.alt_eve.trialinfo = bcc_dat.(bcc_dat.data_name{iE}).trialinfo;
        if ~all(bcc_dat.(bcc_dat.data_name{iE}).trialinfo == cell2mat(sbj_log(:,1))); error([fcfg.sbj_nme ': task list & log do not match!']); end
        
        bcc_dat.(bcc_dat.data_name{iE}).cfg.alt_eve.ovr_all_eve = ones(numel(bcc_dat.(bcc_dat.data_name{iE}).cfg.alt_eve.trialinfo),1)*101;
        
        bcc_dat.(bcc_dat.data_name{iE}).cfg.alt_eve.vis_idn = sbj_lst(sbj_eve_ind,2);
        bcc_dat.(bcc_dat.data_name{iE}).cfg.alt_eve.aud_idn = sbj_lst(sbj_eve_ind,3);
        
        bcc_dat.(bcc_dat.data_name{iE}).cfg.alt_eve.vis_con = cell2mat(sbj_lst(sbj_eve_ind,5));
            bcc_dat.(bcc_dat.data_name{iE}).cfg.alt_eve.vis_con(bcc_dat.(bcc_dat.data_name{iE}).trialinfo==3) = NaN;
        bcc_dat.(bcc_dat.data_name{iE}).cfg.alt_eve.vis_vow = cell2mat(sbj_lst(sbj_eve_ind,6));
            bcc_dat.(bcc_dat.data_name{iE}).cfg.alt_eve.vis_vow(bcc_dat.(bcc_dat.data_name{iE}).trialinfo==3) = NaN;
        
        bcc_dat.(bcc_dat.data_name{iE}).cfg.alt_eve.aud_con = cell2mat(sbj_lst(sbj_eve_ind,7));
            bcc_dat.(bcc_dat.data_name{iE}).cfg.alt_eve.aud_con(bcc_dat.(bcc_dat.data_name{iE}).trialinfo==4) = NaN;
        bcc_dat.(bcc_dat.data_name{iE}).cfg.alt_eve.aud_vow = cell2mat(sbj_lst(sbj_eve_ind,8));
            bcc_dat.(bcc_dat.data_name{iE}).cfg.alt_eve.aud_vow(bcc_dat.(bcc_dat.data_name{iE}).trialinfo==4) = NaN;
        
    else
        
        error(['number of trials don''t match ' cfg.sbj])
        
    end
    
end

cfg = [];
cfg.return_events = 0;
cfg.old_events  = {[1 2] 3};
cfg.new_events  = {101 102};
cfg.crt_alt_eve = 'vis_tot_nse';
cfg.use_alt_eve = 'trialinfo';
bcc_dat = ft_func(@ft_redefine_events,cfg,bcc_dat);

cfg = [];
cfg.return_events = 0;
cfg.old_events  = {[1 2] 4};
cfg.new_events  = {111 112};
cfg.crt_alt_eve = 'aud_tot_nse';
cfg.use_alt_eve = 'trialinfo';
bcc_dat = ft_func(@ft_redefine_events,cfg,bcc_dat);

%% Visual update
eve_hld = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.tsk{1} '_vis.csv'],'[,\t]');

eve_lab = {'UN2_C'   'Orth'    };
cur_lab = {'vis_big' 'vis_ngh' };

for iD = 1:numel(bcc_dat.data_name)
    for iT = 1:numel(sbj_eve_ind)
        
        ort_ind = find(strcmpi(eve_hld(2:end,1),bcc_dat.(bcc_dat.data_name{iE}).cfg.alt_eve.vis_idn{iT})) + 1;
        
        for iE = 1:numel(eve_lab)
            bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.(cur_lab{iE})(iT,1) = eve_hld{ort_ind,strcmpi(eve_hld(1,:),eve_lab{iE})};
        end
        
        if eve_hld{ort_ind,strcmpi(eve_hld(1,:),'FREQ')} > 1; wrd = 101; else wrd = 102; end
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_wrd(iT,1) = wrd;
        
    end
        
end

%
for iD = 1:numel(bcc_dat.data_name)
    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_big(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_wrd==101) = nan;
    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_ngh(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_wrd==101) = nan;
    
    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_big(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo==3) = nan;
    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_ngh(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo==3) = nan;
    
    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_wrd(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo==3) = nan;
end

%
cfg = [];
cfg.eve_rul     = {'med_spl'             'med_spl' };
cfg.new_events  = {[121 121 nan 122 122] [131 131 nan 132 132]};
cfg.crt_alt_eve = {'vis_big_med'         'vis_ngh_med'          };
cfg.fld_nme     = {'vis_big'             'vis_ngh'              };
bcc_dat         = ft_func(@ft_redefine_events2,cfg,bcc_dat);

%% Auditory update
eve_hld = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.tsk{1} '_aud_stt.csv'],'[,\t]');

eve_lab = {'UN2_C'   'Orth'};
cur_lab = {'aud_big' 'aud_ngh'};

for iD = 1:numel(bcc_dat.data_name)
    for iT = 1:numel(sbj_eve_ind)
        
        aud_stm = bcc_dat.(bcc_dat.data_name{iE}).cfg.alt_eve.aud_idn{iT}(3:5);
        if strcmpi(aud_stm(end),'.') || strcmpi(aud_stm(end),'_'); aud_stm = [aud_stm(1:2) 'h']; end
        
        aud_ind = find(strcmpi(eve_hld(2:end,1),aud_stm)) + 1;
        
        for iE = 1:numel(eve_lab)
            bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.(cur_lab{iE})(iT,1) = eve_hld{aud_ind,strcmpi(eve_hld(1,:),eve_lab{iE})};
        end
        
        if eve_hld{aud_ind,strcmpi(eve_hld(1,:),'FREQ')} > 1; wrd = 201; else wrd = 202; end
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_wrd(iT,1) = wrd;
        
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_arp{iT,1} = eve_hld{aud_ind,strcmpi(eve_hld(1,:),'ARPABET')};
        
    end
end

%
for iD = 1:numel(bcc_dat.data_name)
    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_big(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_wrd==201) = nan;
    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_ngh(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_wrd==201) = nan;
    
    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_big(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo==4) = nan;
    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_ngh(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo==4) = nan;
    
    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_wrd(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo==4) = nan;
end

%
cfg = [];
cfg.eve_rul     = {'med_spl'             'med_spl' };
cfg.new_events  = {[221 221 nan 222 222] [231 231 nan 232 232]};
cfg.crt_alt_eve = {'aud_big_med'         'aud_ngh_med'          };
cfg.fld_nme     = {'aud_big'             'aud_ngh'              };
bcc_dat         = ft_func(@ft_redefine_events2,cfg,bcc_dat);

%% Save Stats
if ~exist([fcfg.clr_fld '/' 'events'],'dir'); mkdir([fcfg.clr_fld '/' 'events']); end
if ~exist([fcfg.clr_fld '/' 'events' '/' fcfg.sbj_nme],'dir'); mkdir([fcfg.clr_fld '/' 'events' '/' fcfg.sbj_nme]); end

eve_hld = fieldnames(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve);
eve_hld = sort(eve_hld);
eve_tbl = num2cell([1:numel(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo)]');
for iEH = 1:numel(eve_hld)
    if ~strcmpi(eve_hld{iEH},'orig_trialinfo')
    if isnumeric(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_hld{iEH})(1))
        eve_tbl(:,iEH+1) = num2cell(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_hld{iEH}));
    else
        eve_tbl(:,iEH+1) = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_hld{iEH});
    end
    end
end

cell2csv([fcfg.clr_fld '/' 'events' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '.csv'],[[{''} eve_hld'] ; eve_tbl])

%%
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.sve_app  = 'app_all';
cfg.filename = [outpath '/' sbj '_overall_data.mat'];
ft_func([],cfg,bcc_dat);

end

%% TESTING
% figure()
% hold on
% scatter(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_big(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_wrd==102), ...
%     bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_ngh(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_wrd==102),'r')
% scatter(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_big(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_wrd==101), ...
%     bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_ngh(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.vis_wrd==101),'g')
% 
% figure()
% hold on
% scatter(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_big(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_wrd==202), ...
%     bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_ngh(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_wrd==202),'r')
% scatter(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_big(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_wrd==201), ...
%     bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_ngh(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.aud_wrd==201),'g')
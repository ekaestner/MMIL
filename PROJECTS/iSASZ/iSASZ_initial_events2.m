function iSASZ_initial_events2(fcfg)

%% Initial Load
sbj  = fcfg.sbj_nme;

fprintf([fcfg.sbj_nme ': Starting initial stats work on %s \n'],sbj)

infile  = [fcfg.fle_out_pth '/' 'epoch_data' '/' sbj '_overall_data.mat'];
outpath = [fcfg.fle_out_pth '/' 'epoch_data' '/'];

cfg = [];
cfg.load = 'yes';
cfg.file = [outpath '/' sbj '_overall_data.mat'];
eve_dat  = ft_func([],cfg);

tsk         = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'tsk'); %mmil_readtext([fcfg.clr_fld '/tsk/' subj]);
tsk         = repmat(tsk,1,2);

end_dir     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'end_dir'); %mmil_readtext([fcfg.clr_fld '/end_dir/' subj]);

for iD = 1:numel(eve_dat.data_name);
    
    dat = eve_dat.(eve_dat.data_name{iD});
    
    dat.cfg.alt_eve.trialinfo = dat.cfg.alt_eve.trialinfo;
    
    %% Add in  trial words
    if  strcmpi(tsk{iD},'sa') || isfield(dat.cfg,'task') && strcmpi(dat.tsk{iD},'sa')
        
        aud_wrd = mmil_readtext([fcfg.clr_fld '/' 'task' '/' 'SAi-actual.csv']);
        
        dat.cfg.alt_eve.stm_idn = upper(aud_wrd(:,2));
        if isfield(dat.cfg.alt_eve,'skp'); dat.cfg.alt_eve.stm_idn(dat.cfg.alt_eve.skp) = []; end
        
        num_eve_hld{iD} = 1:numel(aud_wrd(:,2)); if isfield(dat.cfg.alt_eve,'skp'); num_eve_hld{iD}(dat.cfg.alt_eve.skp) = []; end
                
    end
    
    if strcmpi(tsk{iD},'sz') || isfield(dat.cfg,'task') && strcmpi(dat.tsk{iD},'sz')
        
        ort_wrd = mmil_readtext([fcfg.clr_fld '/' 'task' '/' 'SZi-actual.csv']);
        
        dat.cfg.alt_eve.stm_idn = upper(ort_wrd(:,2));
        if isfield(dat.cfg.alt_eve,'skp'); dat.cfg.alt_eve.stm_idn(dat.cfg.alt_eve.skp) = []; end
        
        num_eve_hld{iD} = 1:numel(ort_wrd(:,2)); if isfield(dat.cfg.alt_eve,'skp'); num_eve_hld{iD}(dat.cfg.alt_eve.skp) = []; end
        
    end
    
    
    
    %% Change initial repeated words
    if  strcmpi(tsk{iD},'sa') || isfield(dat.cfg,'task') && strcmpi(dat.tsk{iD},'sa')
        ini_wrd = unique(dat.cfg.alt_eve.stm_idn(ismember(dat.trialinfo,15:18)));
    elseif  strcmpi(tsk{iD},'sz') || isfield(dat.cfg,'task') && strcmpi(dat.tsk{iD},'sz')
        ini_wrd = unique(dat.cfg.alt_eve.stm_idn(ismember(dat.trialinfo,5:8)));
    end
    
    for iRP = 1:numel(ini_wrd); fst_trl(iRP) = find(strcmpi(ini_wrd{iRP},dat.cfg.alt_eve.stm_idn),1); end
    dat.trialinfo(fst_trl) = dat.trialinfo(fst_trl) - 4;
    
    %% Take note of distance between repeated words
    if any(dat.cfg.alt_eve.trialinfo < 10)
        rep_wrd = unique(dat.cfg.alt_eve.stm_idn(dat.cfg.alt_eve.trialinfo==5 | dat.cfg.alt_eve.trialinfo==6  | dat.cfg.alt_eve.trialinfo==7 | dat.cfg.alt_eve.trialinfo==8));
    elseif any(dat.cfg.alt_eve.trialinfo > 10)
        rep_wrd = unique(dat.cfg.alt_eve.stm_idn(dat.cfg.alt_eve.trialinfo==15 | dat.cfg.alt_eve.trialinfo==16  | dat.cfg.alt_eve.trialinfo==17 | dat.cfg.alt_eve.trialinfo==18));
    end
    
    dat.cfg.alt_eve.rep_pos = zeros(numel(dat.cfg.alt_eve.stm_idn),1);
    for iRW = 1:numel(rep_wrd)
        wrd_pos = find(strcmpi(rep_wrd{iRW},dat.cfg.alt_eve.stm_idn));
        dat.cfg.alt_eve.rep_pos(wrd_pos(2:end),1) = diff(wrd_pos);
    end
    
    dat.cfg.alt_eve.rep_pos(dat.cfg.alt_eve.rep_pos==0) = nan;
    
    %% Add in orthographic and phonetic statistics
    if  strcmpi(tsk{iD},'sa') || isfield(dat.cfg,'task') && strcmpi(dat.tsk{iD},'sa')
        phn_stt = mmil_readtext([fcfg.clr_fld '/' 'task' '/' 'sa_phonetic.csv']);
        for iPS = 4:size(phn_stt,2)
            dat.cfg.alt_eve.(phn_stt{1,iPS}) = cell2mat(phn_stt(2:end,iPS));
            if isfield(dat.cfg.alt_eve,'skp'); dat.cfg.alt_eve.(phn_stt{1,iPS})(dat.cfg.alt_eve.skp) = []; end
            dat.cfg.alt_eve.(phn_stt{1,iPS})(dat.cfg.alt_eve.(phn_stt{1,iPS})==-999) = nan;
        end
    end
    
    if strcmpi(tsk{iD},'sz') || isfield(dat.cfg,'task') && strcmpi(dat.tsk{iD},'sz')
        ort_stt = mmil_readtext([fcfg.clr_fld '/' 'task' '/' 'sz_orthographic.csv']);
        for iPS = 4:size(ort_stt,2)
            dat.cfg.alt_eve.(ort_stt{1,iPS}) = cell2mat(ort_stt(2:end,iPS));
            if isfield(dat.cfg.alt_eve,'skp'); dat.cfg.alt_eve.(ort_stt{1,iPS})(dat.cfg.alt_eve.skp) = []; end
            dat.cfg.alt_eve.(ort_stt{1,iPS})(dat.cfg.alt_eve.(ort_stt{1,iPS})==-999) = nan;
        end
    end
    
    %% Add in semantic priming
    if  strcmpi(tsk{iD},'sa') || isfield(dat.cfg,'task') && strcmpi(dat.tsk{iD},'sa')
        aud_sem = mmil_readtext([fcfg.clr_fld '/' 'task' '/' 'aud_prm_val.csv']);
        dat.cfg.alt_eve.sem_prm = cell2mat(aud_sem);
        if isfield(dat.cfg.alt_eve,'skp'); dat.cfg.alt_eve.sem_prm(dat.cfg.alt_eve.skp) = []; end
    end
    
    if strcmpi(tsk{iD},'sz') || isfield(dat.cfg,'task') && strcmpi(dat.tsk{iD},'sz')
        vis_sem = mmil_readtext([fcfg.clr_fld '/' 'task' '/' 'vis_prm_val.csv']);
        dat.cfg.alt_eve.sem_prm = cell2mat(vis_sem);
        if isfield(dat.cfg.alt_eve,'skp'); dat.cfg.alt_eve.sem_prm(dat.cfg.alt_eve.skp) = []; end
    end
    
    %% Choose controlled animal/object
    sem = mmil_readtext([fcfg.clr_fld '/' 'task' '/' 'semantic_database.csv']);
    
    if  strcmpi(tsk{iD},'sa') || isfield(dat.cfg,'task') && strcmpi(dat.tsk{iD},'sa')
        for iPS = [6 8]
            for iT = 1:size(dat.cfg.alt_eve.stm_idn)
                ttt = find(strcmpi(dat.cfg.alt_eve.stm_idn{iT},sem(:,1)));
                if ~isempty(ttt)
                    dat.cfg.alt_eve.(sem{1,iPS})(iT,1) = cell2mat(sem(ttt,iPS));
                else
                    dat.cfg.alt_eve.(sem{1,iPS})(iT,1) = nan;
                end
            end
        end
    end
    
    if strcmpi(tsk{iD},'sz') || isfield(dat.cfg,'task') && strcmpi(dat.tsk{iD},'sz')
        for iPS = [6 8]
            for iT = 1:size(dat.cfg.alt_eve.stm_idn)
                ttt = find(strcmpi(dat.cfg.alt_eve.stm_idn{iT},sem(:,1)));
                if ~isempty(ttt)
                    dat.cfg.alt_eve.(sem{1,iPS})(iT,1) = cell2mat(sem(ttt,iPS));
                else
                    dat.cfg.alt_eve.(sem{1,iPS})(iT,1) = nan;
                end
            end
        end
    end
    
    %% Add in arousal and valence statistics
    emo = mmil_readtext([fcfg.clr_fld '/' 'task' '/' 'EmotionRatings_trimmed.csv']);
    
    if  strcmpi(tsk{iD},'sa') || isfield(dat.cfg,'task') && strcmpi(dat.tsk{iD},'sa')
        for iPS = 3:4
            for iT = 1:size(dat.cfg.alt_eve.stm_idn)
                ttt = find(strcmpi(dat.cfg.alt_eve.stm_idn{iT},emo(:,2)));
                if ~isempty(ttt)
                    dat.cfg.alt_eve.(emo{1,iPS})(iT,1) = cell2mat(emo(ttt,iPS));
                else
                    dat.cfg.alt_eve.(emo{1,iPS})(iT,1) = nan;
                end
            end
        end
    end
    
    if strcmpi(tsk{iD},'sz') || isfield(dat.cfg,'task') && strcmpi(dat.tsk{iD},'sz')
        for iPS = 3:4
            for iT = 1:size(dat.cfg.alt_eve.stm_idn)
                ttt = find(strcmpi(dat.cfg.alt_eve.stm_idn{iT},emo(:,2)));
                if ~isempty(ttt)
                    dat.cfg.alt_eve.(emo{1,iPS})(iT,1) = cell2mat(emo(ttt,iPS));
                else
                    dat.cfg.alt_eve.(emo{1,iPS})(iT,1) = nan;
                end
            end
        end
    end
    
    eve_dat.(eve_dat.data_name{iD}) = dat;
    
end

%% Update Events
if any(strcmpi(tsk,'sz'))
    
    % overall vis
    cfg = [];
    cfg.data_name   = find(strcmpi(tsk,'sz'));
    cfg.eve_rul     = {'new'};
    cfg.old_events  = {{[1 2 3 4 5 6 7 8]}};
    cfg.new_events  = {101};
    cfg.crt_alt_eve = {'ovr'};
    cfg.use_alt_eve = {'trialinfo'};
    eve_dat = ft_func(@ft_redefine_events2,cfg,eve_dat);
    
    % old/new eve % ani/obj
    cfg = [];
    cfg.data_name   = find(strcmpi(tsk,'sz'));
    cfg.eve_rul     = {'new' 'new'};
    cfg.old_events  = {{[1 2 3 4] [5 6 7 8]} {[1 3 5 7] [2 4 6 8]}};
    cfg.new_events  = {[111 112]              [121 122]};
    cfg.crt_alt_eve = {'new_old' 'ani_obj'};
    cfg.use_alt_eve = {'trialinfo'   'trialinfo'};
    eve_dat = ft_func(@ft_redefine_events2,cfg,eve_dat);
    
    % old/new close
    cfg = [];
    cfg.data_name   = find(strcmpi(tsk,'sz'));
    cfg.eve_rul     = {{'all' '<=10'}};
    cfg.fld_nme     = {'rep_pos'};
    cfg.old_events  = {{[1 2 3 4] [5 6 7 8]}};
    cfg.new_events  = {[113 114]};
    cfg.crt_alt_eve = {'new_old_sht'};
    cfg.use_alt_eve = {'trialinfo' };
    eve_dat = ft_func(@ft_redefine_events2,cfg,eve_dat);
    
    % ani/obj priming    
    cfg = [];
    cfg.data_name   = find(strcmpi(tsk,'sz'));
    cfg.eve_rul     = {'med_spl'};
    cfg.new_events  = {[1011 1011 0 0 1012 1012]};
    cfg.crt_alt_eve = {'cat_prm_new'};
    cfg.fld_nme = {'sem_prm'};
    eve_dat = ft_func(@ft_redefine_events2,cfg,eve_dat);
    
    % density % freq % Bigram % Length
    cfg = [];
    cfg.data_name   = find(strcmpi(tsk,'sz'));
    cfg.eve_rul     = {'med_spl' 'med_spl' 'med_spl'};
    cfg.new_events  = {[131 131 0 132 132] [141 141 0 142 142] [151 151 0 152 152]};
    cfg.crt_alt_eve = {'LEN_med' 'UN2_C_med' 'FREQ_med'};
    cfg.fld_nme     = {'LEN' 'UN2_C' 'FREQ'};
    eve_dat = ft_func(@ft_redefine_events2,cfg,eve_dat);
    
    % Arousal % Valence
    cfg = [];
    cfg.data_name   = find(strcmpi(tsk,'sz'));
    cfg.eve_rul     = {'med_spl'           'med_spl'           'cmb'};
    cfg.new_events  = {[161 161 nan 162 162] [171 171 nan 172 172] [181 182 183 184]};
    cfg.old_events  = {[]                  []                  {[161 171] [161 172] [162 171] [162 172]}};
    cfg.crt_alt_eve = {'val'               'ars'               'emo'};
    cfg.use_alt_eve = {''                  ''                  {'val' 'ars'}};
    cfg.fld_nme     = {'Val_Mean_Sum'      'Ars_Mean_Sum'      ''};
    eve_dat = ft_func(@ft_redefine_events2,cfg,eve_dat);
    
    % Semantics
    cfg = [];
    cfg.data_name   = find(strcmpi(tsk,'sz'));
    cfg.eve_rul     = {'med_spl'           'med_spl'           'cmb'};
    cfg.new_events  = {[191 191 nan 192 192] [193 193 nan 194 194] [195 196 197 198]};
    cfg.old_events  = {[]                  []                  {[191 193] [191 194] [192 193] [192 194]}};
    cfg.crt_alt_eve = {'sem_mnp'      'sem_mot'                'sem_spl'};
    cfg.use_alt_eve = {''             ''                       {'sem_mnp' 'sem_mot'}};
    cfg.fld_nme     = {'manipulation' 'motion'      ''};
    eve_dat = ft_func(@ft_redefine_events2,cfg,eve_dat);
    
end

if any(strcmpi(tsk,'sa'))
    
    % overall aud
    cfg = [];
    cfg.data_name   = find(strcmpi(tsk,'sa'));
    cfg.eve_rul     = {'new'};
    cfg.old_events  = {{[11 12 13 14 15 16 17 18]}};
    cfg.new_events  = {201};
    cfg.crt_alt_eve = {'ovr'};
    cfg.use_alt_eve = {'trialinfo'};
    eve_dat = ft_func(@ft_redefine_events2,cfg,eve_dat);
    
    % old/new eve % ani/obj
    cfg = [];
    cfg.data_name   = find(strcmpi(tsk,'sa'));
    cfg.eve_rul     = {'new' 'new'};
    cfg.old_events  = {{[11 12 13 14] [15 16 17 18]} {[11 13 15 17] [12 14 16 18]}};
    cfg.new_events  = {[211 212]                     [221 222]};
    cfg.crt_alt_eve = {'new_old' 'ani_obj'};
    cfg.use_alt_eve = {'trialinfo'   'trialinfo'};
    eve_dat = ft_func(@ft_redefine_events2,cfg,eve_dat);
       
    % old/new close
    cfg = [];
    cfg.data_name   = find(strcmpi(tsk,'sa'));
    cfg.eve_rul     = {{'all' '<= 15'}};
    cfg.fld_nme     = {'rep_pos'};
    cfg.old_events  = {{[11 12 13 14] [15 16 17 18]}};
    cfg.new_events  = {[213 214]};
    cfg.crt_alt_eve = {'new_old_sht'};
    cfg.use_alt_eve = {'trialinfo' };
    eve_dat = ft_func(@ft_redefine_events2,cfg,eve_dat);
    
    % ani/obj priming
    cfg = [];
    cfg.data_name   = find(strcmpi(tsk,'sa'));
    cfg.eve_rul     = {'med_spl'};
    cfg.new_events  = {[2011 2011 0 0 2012 2012]};
    cfg.crt_alt_eve = {'cat_prm_new'};
    cfg.fld_nme = {'sem_prm'};
    eve_dat = ft_func(@ft_redefine_events2,cfg,eve_dat);
    
    % density % freq % Bigram % Length
    cfg = [];
    cfg.data_name   = find(strcmpi(tsk,'sa'));
    cfg.eve_rul     = {'med_spl' 'med_spl' 'med_spl'};
    cfg.new_events  = {[231 231 nan 232 232] [241 241 nan 242 242] [251 251 nan 252 252]};
    cfg.crt_alt_eve = {'LEN_med'           'UN2_C_med'         'FREQ_med'};
    cfg.fld_nme     = {'LEN'               'UN2_C'             'FREQ'    };
    eve_dat = ft_func(@ft_redefine_events2,cfg,eve_dat);
    
    % Arousal % Valence
    cfg = [];
    cfg.data_name   = find(strcmpi(tsk,'sa'));
    cfg.eve_rul     = {'med_spl' 'med_spl'};
    cfg.new_events  = {[271 271 nan 272 272] [281 281 nan 282 282]};
    cfg.crt_alt_eve = {'Valence'           'Arousal'};
    cfg.fld_nme     = {'Val_Mean_Sum'      'Ars_Mean_Sum'};
    eve_dat = ft_func(@ft_redefine_events2,cfg,eve_dat);
    
    % Arousal % Valence
    cfg = [];
    cfg.data_name   = find(strcmpi(tsk,'sa'));
    cfg.eve_rul     = {'med_spl'           'med_spl'           'cmb'};
    cfg.new_events  = {[261 261 nan 262 262] [271 271 nan 272 272] [281 282 283 284]};
    cfg.old_events  = {[]                  []                  {[261 271] [261 272] [262 271] [262 272]}};
    cfg.crt_alt_eve = {'val'               'ars'           'emo'};
    cfg.use_alt_eve = {''                  ''                  {'val' 'ars'}};
    cfg.fld_nme     = {'Val_Mean_Sum'      'Ars_Mean_Sum'      ''};
    eve_dat = ft_func(@ft_redefine_events2,cfg,eve_dat);
    
    % Semantics
    cfg = [];
    cfg.data_name   = find(strcmpi(tsk,'sa'));
    cfg.eve_rul     = {'med_spl'           'med_spl'           'cmb'};
    cfg.new_events  = {[291 291 nan 292 292] [293 293 nan 294 294] [295 296 297 298]};
    cfg.old_events  = {[]                  []                  {[291 293] [291 294] [292 293] [292 294]}};
    cfg.crt_alt_eve = {'sem_mnp'      'sem_mot'                'sem_spl'};
    cfg.use_alt_eve = {''             ''                       {'sem_mnp' 'sem_mot'}};
    cfg.fld_nme     = {'manipulation' 'motion'      ''};
    eve_dat = ft_func(@ft_redefine_events2,cfg,eve_dat);
    
end

%% Save Stat Output
if ~exist([fcfg.clr_fld '/' 'events' '/' fcfg.sbj_nme]); mkdir([fcfg.clr_fld '/' 'events' '/' fcfg.sbj_nme]); end

for iD = 1:numel(eve_dat.data_name);
    stt_hld = num2cell(num_eve_hld{iD}');
    stt_nme = fieldnames(eve_dat.(eve_dat.data_name{iD}).cfg.alt_eve);
    stt_nme(strcmpi(stt_nme,'skp')) = [];
    for iST = 1:numel(stt_nme)
        if all(isnumeric(eve_dat.(eve_dat.data_name{iD}).cfg.alt_eve.(stt_nme{iST})));
        stt_hld(:,iST+1) = num2cell(eve_dat.(eve_dat.data_name{iD}).cfg.alt_eve.(stt_nme{iST}));
        else
            stt_hld(:,iST+1) = eve_dat.(eve_dat.data_name{iD}).cfg.alt_eve.(stt_nme{iST});
        end
    end  
    cell2csv([fcfg.clr_fld '/' 'events' '/' fcfg.sbj_nme '/' eve_dat.data_name{iD}],[['' '' fieldnames(eve_dat.(eve_dat.data_name{iD}).cfg.alt_eve)]';stt_hld]);
    
end

%% Save Data & Stats
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.sve_app  = 'app_all';
cfg.filename =[outpath '/' sbj '_overall_data.mat'];
ft_func([],cfg,eve_dat);

end
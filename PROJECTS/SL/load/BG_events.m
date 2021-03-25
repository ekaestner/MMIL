function BG_events(fcfg)

%% Load
cfg = [];
cfg.load = 'yes';
cfg.file = [fcfg.out_pth '/' fcfg.sbj_nme '_overall_data.mat'];
bgr_dat  = ft_func([],cfg);

%%
tsk = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'tsk'); %mmil_readtext([fcfg.clr_fld '/tsk/' subj]);     

lst = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_' tsk{1} '.csv']);

end_ind = find(all(cellfun(@isempty,lst)));

nme = {'stm_idn' 'trialinfo' 'lng_nse'};
ind = [1          2           3];

for iD = 1:numel(bgr_dat.data_name)
    for iE = 1:numel(nme)
        
        if isnumeric(lst{1,ind(iE)})
            bgr_dat.(bgr_dat.data_name{iD}).cfg.alt_eve.(nme{iE}) = cell2mat(lst(:,ind(iE)));
        else
            bgr_dat.(bgr_dat.data_name{iD}).cfg.alt_eve.(nme{iE}) = lst(:,ind(iE));
        end
        
        if numel(bgr_dat.(bgr_dat.data_name{iD}).trialinfo) ~= numel(bgr_dat.(bgr_dat.data_name{iD}).cfg.alt_eve.(nme{iE})); error('ERROR ERROR EXTERMINATE EXTERMINATE'); end
        
    end
    
    vis_ind = bgr_dat.(bgr_dat.data_name{iD}).cfg.alt_eve.lng_nse==1  |  bgr_dat.(bgr_dat.data_name{iD}).cfg.alt_eve.lng_nse==2;
    aud_ind = bgr_dat.(bgr_dat.data_name{iD}).cfg.alt_eve.lng_nse==11 |  bgr_dat.(bgr_dat.data_name{iD}).cfg.alt_eve.lng_nse==12;
    
    bgr_dat.(bgr_dat.data_name{iD}).cfg.alt_eve.ovr(vis_ind,1) = 101;
    bgr_dat.(bgr_dat.data_name{iD}).cfg.alt_eve.ovr(aud_ind,1) = 201;
    
end

%% Save
cfg = [];
cfg.str_nme  = 'bgr_dat';
cfg.save     = 'yes';
cfg.filename =[fcfg.out_pth '/' fcfg.sbj_nme '_overall_data.mat'];
ft_func([],cfg,bgr_dat);

end
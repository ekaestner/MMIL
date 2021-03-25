%% Count Channels
clear; clc

outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/paper_plot/ThreeStat/sig_chn';
inpath  = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/paper_plot/ThreeStat/sig_chn';
sbj_cnt = 1;

% Setting up Variables
sbj_hld  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/subjects');

dat_typ = {'lfp' 'hgp' 'tht' 'bta'};
ele_hld = {'subject','v_act_ely','v_rep_ely','a_act_ely', 'a_rep_ely','ovr_lap_ovr','ovr_lap_rep'};

for sbj_num = [1:17 19:numel(sbj_hld)-1]; 
    
    subj  = sbj_hld{sbj_num};
    
    cfg = [];
    cfg.load = 'yes';
    cfg.file = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/' subj '_overall_data.mat'];
    sem_dat  = ft_func([],cfg);
    
    sbj_cnt = sbj_cnt + 1;
    
    has_vis = any(sem_dat.(sem_dat.data_name{1}).trialinfo < 9);
    has_aud = any(sem_dat.(sem_dat.data_name{1}).trialinfo > 9);
    
    for iMV = 1:numel(sem_dat.data_name)
        
        % Make list of overall numbers of electrodes
        sig_chn     = mmil_readtext([inpath '/' subj '/' sem_dat.data_name{iMV}]);
        sig_chn_int = sig_chn(2:end,[1 2 4 8 10]);
        
        ttt.(dat_typ{iMV}){sbj_cnt,1} = subj;
        
        ind_icd = 1:size(sig_chn_int,1);
        dpt_ind = find(~cellfun(@isempty,strfind(sig_chn_int(:,1),'D')) | ~cellfun(@isempty,strfind(sig_chn_int(:,1),'d')));
        if any(~cellfun(@isempty,strfind(sig_chn_int(dpt_ind,1),'>')))
            str_ind = strfind(sig_chn_int(dpt_ind,1),'>');
        else
            str_ind = strfind(sig_chn_int(dpt_ind,1),'-');
        end
        if ~isempty(dpt_ind); dpt_ind = dpt_ind(cellfun(@(x,y) strcmpi(x(y+1),'d'),sig_chn_int(dpt_ind,1),str_ind)); end
        ind_icd(dpt_ind-1) = [];
        
        if has_vis && has_aud
            
            ttt.(dat_typ{iMV}){sbj_cnt,2} = sum(cell2mat(sig_chn_int(ind_icd,2)));
            ttt.(dat_typ{iMV}){sbj_cnt,3} = sum(cell2mat(sig_chn_int(ind_icd,3)));
            ttt.(dat_typ{iMV}){sbj_cnt,4} = sum(cell2mat(sig_chn_int(ind_icd,4)));
            ttt.(dat_typ{iMV}){sbj_cnt,5} = sum(cell2mat(sig_chn_int(ind_icd,5)));
            ttt.(dat_typ{iMV}){sbj_cnt,6} = sum(cell2mat(sig_chn_int(ind_icd,4)) & cell2mat(sig_chn_int(ind_icd,2)));
            ttt.(dat_typ{iMV}){sbj_cnt,7} = sum(cell2mat(sig_chn_int(ind_icd,5)) & cell2mat(sig_chn_int(ind_icd,3)));
            
        elseif has_vis
            
            ttt.(dat_typ{iMV}){sbj_cnt,2} = sum(cell2mat(sig_chn_int(ind_icd,2)));
            ttt.(dat_typ{iMV}){sbj_cnt,3} = sum(cell2mat(sig_chn_int(ind_icd,3)));
            
        elseif has_aud
            
            ttt.(dat_typ{iMV}){sbj_cnt,4} = sum(cell2mat(sig_chn_int(ind_icd,2)));
            ttt.(dat_typ{iMV}){sbj_cnt,5} = sum(cell2mat(sig_chn_int(ind_icd,3)));
            
        end
        
    end
    
end

for iMV = 1:numel(dat_typ)
   
    ele_hld_typ = [ele_hld ; ttt.(dat_typ{iMV})];
    cell2csv(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/paper_plot/ThreeStat/sig_chn/' dat_typ{iMV} '.csv'],ele_hld_typ);
    
end
dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data';

sbj_nme = 'NY007_SA_SZ';

% Load Data
cfg = [];
cfg.load = 'yes';
cfg.file = [dat_fld '/' sbj_nme '_overall_data.mat'];
ecg_dat  = ft_func([],cfg);

% Get trial count
vis_new = find(ecg_dat.(ecg_dat.data_name{2}).cfg.alt_eve.vis_new_old==111); vis_new_num = numel(vis_new);
vis_old = find(ecg_dat.(ecg_dat.data_name{2}).cfg.alt_eve.vis_new_old==112); vis_old_num = numel(vis_old);

aud_new = find(ecg_dat.(ecg_dat.data_name{2}).cfg.alt_eve.aud_new_old==211); aud_new_num = numel(aud_new);
aud_old = find(ecg_dat.(ecg_dat.data_name{2}).cfg.alt_eve.aud_new_old==212); aud_old_num = numel(aud_new);

if aud_new > vis_new && aud_old > vis_old
    
    lng_fin = 1; ngh_fin = 1; frq_fin = 1;
    
    while lng_fin || ngh_fin || frq_fin
        
        % New Numbers %%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Get Wrd Frq, Orth Neigh, Length, Phoneme Length
        vis_lng = ecg_dat.(ecg_dat.data_name{2}).cfg.alt_eve.vis_lng(vis_new);
        vis_frq = ecg_dat.(ecg_dat.data_name{2}).cfg.alt_eve.vis_frq(vis_new);
        vis_ngh = ecg_dat.(ecg_dat.data_name{2}).cfg.alt_eve.vis_ngh(vis_new);
        
        aud_chs_new = datasample(aud_new,vis_new_num); % Choose randomly for equal numbers
        
        aud_lng = ecg_dat.(ecg_dat.data_name{2}).cfg.alt_eve.aud_lng(aud_chs_new);
        aud_frq = ecg_dat.(ecg_dat.data_name{2}).cfg.alt_eve.aud_frq(aud_chs_new);
        aud_ngh = ecg_dat.(ecg_dat.data_name{2}).cfg.alt_eve.aud_ngh(aud_chs_new);     
        
        % ttest Wrd Frq, Orth Neigh, Length, Phoneme Length
        [~,frq_pvl] = ttest(vis_frq,aud_frq); %[frq_pvl nanmean(vis_frq) nanmean(aud_frq)]
            if frq_pvl>.20; frq_fin = 0; else frq_fin = 1; end
        [~,ngh_pvl] = ttest(vis_ngh,aud_ngh); %[ngh_pvl nanmean(vis_ngh) nanmean(aud_ngh)]
            if ngh_pvl>.20; ngh_fin = 0; else ngh_fin = 1; end
        [~,lng_pvl] = ttest(vis_lng,aud_lng); %[lng_pvl nanmean(vis_lng) nanmean(aud_lng)]
            if lng_pvl>.20; lng_fin = 0; else lng_fin = 1; end
        
        % Old Numbers %%%%%%%%%%%%%%%%%%%%%%%%%%%
          
        
        
        
        
    end
        
end
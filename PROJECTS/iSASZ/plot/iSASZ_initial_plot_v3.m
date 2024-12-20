function iSASZ_initial_plot_v3(fcfg)

fprintf(['Plotting: ' fcfg.sbj_nme])

%% Setup
cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.fle_out_pth '/' 'epoch_data' '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

for iD = 1:numel(bcc_dat.data_name); bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_lab.stt_lab = bcc_dat.(bcc_dat.data_name{iD}).label; bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_lab.label = bcc_dat.(bcc_dat.data_name{iD}).label; end

if any(bcc_dat.(bcc_dat.data_name{1}).trialinfo<10)
    
    dat_loc(:,1:2) = ones(4,2);
    
    alt_eve(:,1:2) = {'vis_ovr'         'vis_lng_med' ; ...
                      'vis_new_old'     'vis_big_med' ; ...
                      'vis_ani_obj'     'vis_frq_med' ; ...
                      'vis_sem_prm_med' 'vis_ngh_med' };
                  
    eve(:,1:2)     = {[101]     [131 132] ; ...
                      [111 112] [141 142] ; ...
                      [121 122] [151 152] ; ...
                      [123 124] [161 162] };
                  
    col_ord(:,1:2) = {{rgb('black')}                            {rgb('dark red') rgb('bright red')} ; ...
                      {rgb('red') rgb('orange')}                {rgb('dark red') rgb('bright red')} ; ...
                      {rgb('green') rgb('grey')}                {rgb('dark red') rgb('bright red')} ; ...
                      {rgb('dark purple') rgb('bright purple')} {rgb('dark red') rgb('bright red')} } ;
                  
    plt_lbl(:,1:2) = {'Visual Overall'          'Visual Length' ; ...
                      'Visual Old / New'        'Visual Bigram' ; ...
                      'Visual Animal / Object'  'Visual Frequency' ; ...
                      'Visual Semantic Priming' 'Visual Neighborhood'};
    
    stt_dat(:,1:2) = {{'vis_ovr_stt'}     {'vis_lng_stt'} ; ...
                      {'vis_new_old_stt'} {'vis_big_stt'} ; ...
                      {'vis_ani_obj_stt'} {'vis_frq_stt'} ; ...
                      {'vis_sem_prm_stt'} {'vis_ngh_stt'} };
                  
    stt_col(:,1:2) = {{rgb('grey')}   {rgb('red')} ; ...
                      {rgb('light orange')} {rgb('red')} ; ...
                      {rgb('green')}  {rgb('red')} ; ...
                      {rgb('purple')} {rgb('red')} };
   
    v_lne(:,1:2)     = {0 0 ; ...
                        0 0 ; ...
                        0 0 ; ...
                        0 0 };
                    
    v_lne_col(:,1:2) = {{rgb('red')} {rgb('red')} ; ...
                        {rgb('red')} {rgb('red')} ; ...
                        {rgb('red')} {rgb('red')} ; ...
                        {rgb('red')} {rgb('red')} };
                  
else
    
    dat_loc(:,1:2) = zeros(4,2);
    alt_eve(:,1:2) = {'' '' ; ...
                      '' '' ; ...
                      '' '' ; ...
                      '' '' };
    eve(:,1:2)     = {[] [] ; ...
                      [] [] ; ...
                      [] [] ; ...
                      [] [] };
    
    col_ord(:,1:2) = {{} {} ;
                      {} {} ; 
                      {} {} ;
                      {} {} };
                  
    plt_lbl(:,1:2) = {'' '' ; ...
                      '' '' ; ...
                      '' '' ; ...
                      '' '' };
    
    stt_dat(:,1:2) = {{} {} ; ...
                      {} {} ; ...
                      {} {} ; ...
                      {} {} };
    stt_col(:,1:2) = {{} {} ; ...
                      {} {} ; ...
                      {} {} ; ...
                      {} {} };
                  
    v_lne(:,1:2)    = {[] [] ; ...
                       [] [] ; ...
                       [] [] ; ...
                       [] [] };
    v_lne_col(:,1:2) = {{} {} ; ...
                        {} {} ; ...
                        {} {} ; ...
                        {} {} }; 
                    
end

if any(bcc_dat.(bcc_dat.data_name{1}).trialinfo>10)
     
    dat_loc(:,3:4) = ones(4,2);
    
    alt_eve(:,3:4) = {'aud_ovr'         'aud_lng_med' ; ...
                      'aud_new_old'     'aud_big_med' ; ...
                      'aud_ani_obj'     'aud_frq_med' ; ...
                      'aud_sem_prm_med' 'aud_ngh_med' };
                  
    eve(:,3:4)     = {[201]     [231 232] ; ...
                      [211 212] [241 242] ; ...
                      [221 222] [251 252] ; ...
                      [223 224] [261 262] };
                  
    col_ord(:,3:4) = {{rgb('black')}                            {rgb('dark blue') rgb('bright blue')} ; ...
                      {rgb('blue') rgb('dark orange')}          {rgb('dark blue') rgb('bright blue')} ; ...
                      {rgb('green') rgb('grey')}                {rgb('dark blue') rgb('bright blue')} ; ...
                      {rgb('dark purple') rgb('bright purple')} {rgb('dark blue') rgb('bright blue')} } ;
                  
    plt_lbl(:,3:4) = {'Auditory Overall'          'Auditory Length' ; ...
                      'Auditory Old / New'        'Auditory Bigram' ; ...
                      'Auditory Animal / Object'  'Auditory Frequency' ; ...
                      'Auditory Semantic Priming' 'Auditory Neighborhood'};
    
    stt_dat(:,3:4) = {{'aud_ovr_stt'}     {'aud_lng_stt'} ; ...
                      {'aud_new_old_stt'} {'aud_big_stt'} ; ...
                      {'aud_ani_obj_stt'} {'aud_frq_stt'} ; ...
                      {'aud_sem_prm_stt'} {'aud_ngh_stt'} };
                  
    stt_col(:,3:4) = {{rgb('grey')}   {rgb('blue')} ; ...
                      {rgb('light orange')} {rgb('blue')} ; ...
                      {rgb('green')}  {rgb('blue')} ; ...
                      {rgb('purple')} {rgb('blue')} };
   
    v_lne(:,3:4)     = {0 0 ; ...
                        0 0 ; ...
                        0 0 ; ...
                        0 0 };
                    
    v_lne_col(:,3:4) = {{rgb('blue')} {rgb('blue')} ; ...
                        {rgb('blue')} {rgb('blue')} ; ...
                        {rgb('blue')} {rgb('blue')} ; ...
                        {rgb('blue')} {rgb('blue')} };
                  
else
    
    dat_loc(:,3:4) = zeros(4,2);
    alt_eve(:,3:4) = {'' '' ; ...
                      '' '' ; ...
                      '' '' ; ...
                      '' '' };
    eve(:,3:4)     = {[] [] ; ...
                      [] [] ; ...
                      [] [] ; ...
                      [] [] };
    
    col_ord(:,3:4) = {{} {} ;
                      {} {} ; 
                      {} {} ;
                      {} {} };
                  
    plt_lbl(:,3:4) = {'' '' ; ...
                      '' '' ; ...
                      '' '' ; ...
                      '' '' };
    
    stt_dat(:,3:4) = {{} {} ; ...
                      {} {} ; ...
                      {} {} ; ...
                      {} {} };
    stt_col(:,3:4) = {{} {} ; ...
                      {} {} ; ...
                      {} {} ; ...
                      {} {} };
                  
    v_lne(:,3:4)    = {[] [] ; ...
                       [] [] ; ...
                       [] [] ; ...
                       [] [] };
    v_lne_col(:,3:4) = {{} {} ; ...
                        {} {} ; ...
                        {} {} ; ...
                        {} {} };            
end

%% Make plots
for iD = 1:2
    
    cfg = [];
    
    cfg.type      = 'chan';
    
    cfg.dat       = {bcc_dat.(bcc_dat.data_name{iD})};
    
    cfg.plt_dim   = [4 4];
    
    cfg.dat_loc   = dat_loc;
    cfg.alt_eve   = alt_eve;
    cfg.eve       = eve;
    
    cfg.plt_lbl   = plt_lbl;
    cfg.ttl_num   = zeros(4,4);
    
    cfg.y_lim     = 'maxmin';
    
    cfg.lgd       = 0;
    cfg.std_err   = 1;
    cfg.lnstyle.col_ord = col_ord;
    cfg.cnd_nme         = {};
    
    cfg.alt_lbl         = 'label';
    
    cfg.stt_dat = stt_dat;
    cfg.stt_col = stt_col;
    cfg.stt_lab = 'stt_lab';
    
    cfg.x_lim       = [-0.3 1.5];
    cfg.v_lne       = v_lne;
    cfg.v_lne_col   = v_lne_col;
    cfg.v_lne_wdt   = 2;
    cfg.axe_fnt_sze = [15 15 15 15 ; ...
        15 15 15 15 ; ...
        15 15 15 15 ; ...
        15 15 15 15 ];
    cfg.axe_lne_sze = [3 3 3 3 ; ...
        3 3 3 3 ; ...
        3 3 3 3 ; ...
        3 3 3 3 ];
    cfg.ttl_lne_sze = [36 36 36 36 ; ...
        36 36 36 36 ; ...
        36 36 36 36 ; ...
        36 36 36 36 ];
    
    cfg.print      = 1;
    cfg.nofig      = 1;
    cfg.print_type = 'png';
    cfg.outdir     = [fcfg.fle_out_pth '/' 'epoch_data' '/' 'initial_plot' '/' fcfg.sbj_nme '/' bcc_dat.data_name{iD}(end-2:end)];
    cfg.prefix     = [fcfg.sbj_nme '_'];
    
    mmil_ieeg_sensor_plot_v5(cfg)
    
end

end
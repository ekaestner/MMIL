clear; clc;

%%
prj_dir = '/home/ekaestne/PROJECTS/';
red_fle = 'mmil_redcap_2019_10_09.csv';

%%
% Redcap Data %%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.prj_dir = prj_dir;
fcfg.red_fle = red_fle;
[ ~ , ~ , ~ , sbj_cog, sbj_emo] = mmil_load_redcap(fcfg);

sbj_nme = sbj_emo.sbj_nme(~isnan(sbj_emo.bai));

fcfg = [];
fcfg.sbj_nme = sbj_nme;
fcfg.prj_dir = prj_dir;
fcfg.red_fle = red_fle;
[sbj_dem , sbj_sze , sbj_scn , sbj_cog, sbj_emo] = mmil_load_redcap(fcfg);

%%
% Total Number
num_con = numel( string_find(sbj_emo.sbj_nme,{'fc'}) );
    ind_con = string_find(sbj_emo.sbj_nme,{'fc'});
num_epd = numel( string_find(sbj_emo.sbj_nme,{'epd'}) );
    ind_epd = string_find(sbj_emo.sbj_nme,{'epd'});

% Correlation of BDI & BAI
subplot(1,2,1)
scatter( sbj_emo.bdi(ind_con) , sbj_emo.bai(ind_con) );
    ylabel('Beck''s Anxiety');xlabel('Beck''s Depression'); xlim([0 63]); ylim([0 63]);
subplot(1,2,2)
scatter( sbj_emo.bdi(ind_epd) , sbj_emo.bai(ind_epd) );
    ylabel('Beck''s Anxiety');xlabel('Beck''s Depression'); xlim([0 63]); ylim([0 63]);
    
% Characterization of BAI
subplot(1,2,1)
hist(sbj_emo.bai(ind_con),50); xlim([0 63]);
subplot(1,2,2)
hist(sbj_emo.bai(ind_epd),50); xlim([0 63]);

con_typ(1) = sum( sbj_emo.bai(ind_con)<=9 );
con_typ(2) = sum( sbj_emo.bai(ind_con)>9 & sbj_emo.bai(ind_con)<=16 );
con_typ(3) = sum( sbj_emo.bai(ind_con)>16 & sbj_emo.bai(ind_con)<=29 );
con_typ(4) = sum( sbj_emo.bai(ind_con)>29 );

epd_typ(1) = sum( sbj_emo.bai(ind_epd)<=9 );
epd_typ(2) = sum( sbj_emo.bai(ind_epd)>9 & sbj_emo.bai(ind_epd)<=16 );
epd_typ(3) = sum( sbj_emo.bai(ind_epd)>16 & sbj_emo.bai(ind_epd)<=29 );
epd_typ(4) = sum( sbj_emo.bai(ind_epd)>29 );
      

% Correlation of BAI with Memory
subplot(1,2,1)
scatter( sbj_cog.log_mem_nor_scr_one(ind_con) , sbj_emo.bai(ind_con) );
    ylabel('Beck''s Anxiety');xlabel('Log Mem 1'); xlim([0 27]); ylim([0 63]);
subplot(1,2,2)
scatter( sbj_cog.log_mem_nor_scr_one(ind_epd) , sbj_emo.bai(ind_epd) );
    ylabel('Beck''s Anxiety');xlabel('Log Mem 1'); xlim([0 27]); ylim([0 63]);

%
subplot(1,2,1)
scatter( sbj_cog.log_mem_nor_scr_two(ind_con) , sbj_emo.bai(ind_con) );
    ylabel('Beck''s Anxiety');xlabel('Log Mem 2'); xlim([0 27]); ylim([0 63]);
subplot(1,2,2)
scatter( sbj_cog.log_mem_nor_scr_two(ind_epd) , sbj_emo.bai(ind_epd) );
    ylabel('Beck''s Anxiety');xlabel('Log Mem 2'); xlim([0 27]); ylim([0 63]);

%
subplot(1,2,1)
scatter( sbj_cog.cvl_lfr_nor_scr(ind_con) , sbj_emo.bai(ind_con) );
    ylabel('Beck''s Anxiety');xlabel('CVLT LFR'); xlim([0 8]); ylim([0 63]);
subplot(1,2,2)
scatter( sbj_cog.cvl_lfr_nor_scr(ind_epd) , sbj_emo.bai(ind_epd) );
    ylabel('Beck''s Anxiety');xlabel('CVLT LFR'); xlim([0 8]); ylim([0 63]);

%
subplot(1,2,1)
scatter( sbj_cog.vp1_nor_scr_pst(ind_con) , sbj_emo.bai(ind_con) );
    ylabel('Beck''s Anxiety');xlabel('VP1'); xlim([0 15]); ylim([0 63]);
subplot(1,2,2)
scatter( sbj_cog.vp1_nor_scr_pst(ind_epd) , sbj_emo.bai(ind_epd) );
    ylabel('Beck''s Anxiety');xlabel('VP1'); xlim([0 15]); ylim([0 63]);

%
subplot(1,2,1)
scatter( sbj_cog.vp2_nor_scr_pst(ind_con) , sbj_emo.bai(ind_con) );
    ylabel('Beck''s Anxiety');xlabel('VP2'); xlim([0 15]); ylim([0 63]);
subplot(1,2,2)
scatter( sbj_cog.vp2_nor_scr_pst(ind_epd) , sbj_emo.bai(ind_epd) );
    ylabel('Beck''s Anxiety');xlabel('VP2'); xlim([0 15]); ylim([0 63]);









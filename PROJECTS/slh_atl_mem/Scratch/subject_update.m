clear; clc;

dta_dir = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/Data/subject_check';
out_dir = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/Data/subject_check/';

%% Emory Files
% Dan
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'Dan_Emory_SLAH_Memory_11.22.21_cleanNONAMES.csv'];
fcfg.dta_col = 2;
[ dan_dta, dan_dta_sbj, dan_dta_col] = ejk_dta_frm( fcfg );

% Rebecca
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'Rebecca_forErik_1_5_23.csv'];
fcfg.dta_col = 2;
[ reb_dta, reb_dta_sbj, reb_dta_col] = ejk_dta_frm( fcfg );
    reb_dta_sbj(cellfun(@isempty,reb_dta_sbj)) = '';   
    
%% Modify Dan Data
clear col_tgh

dan_dta(:,end+1) = dan_dta(:,strcmpi(dan_dta_col,'Surgery'));
dan_dta_col      = [ dan_dta_col {'SurgerySide'} ];

col_tgh{find(strcmpi(dan_dta_col,'Surgery'))} = { 'ATL'  'Right ATL'  ; ...
                                                  'ATL'  'Right ATL (limited open)'  ; ...
                                                  ...
                                                  'SLAH' 'Left SLAH' ; ...
                                                  'SLAH' 'Right SLAH' ; ...
                                                  'SLAH' 'left SLAH' ; ...
                                                  ...
                                                  'Multiple' 'Left SLAH (cleaned up left SLAH after 6 months)' ; ...
                                                  'Multiple' 'Left SLAH with repeat SLAH for more complete ablation' ; ...
                                                  'Multiple' 'Left SLAH- pt. suffered a surgical complication involving a bleed' ; ...
                                                  'Multiple' 'Two right SLAH procedures (memory declined after 2nd - took broader medial TL structures' ; ...
                                                  'Multiple' 'right SLAH; right ATL; right FL resection' ; ...
                                                  ...
                                                  'Other' 'Left RNS' ; ...
                                                  'Other' 'Left SAH (open)' ; ...
                                                  'Other' 'Left SLA involving encephalocele in left temporal lobe' ; ...
                                                  'Other' 'Right TL open resection (not standard)' ; ...
                                                  'Other' 'left SAH (open)' ; ...
                                                  'Other' 'left Temporal lobe RNS' ; ...
                                                  'Other' 'right SLAH (temporal pole extending amygdla; sparring hippocampus)' ; ...
                                                  'Other' 'right superiror temporal gyrus RF ablation' ; ...
                                                  ...
                                                  'None' 'empty' ; ...
                                                  'None' 'no surgery' ; ...
                                                  'None' 'no surgery up to date' ; ...
                                                  'None' 'no surgery-just had SEEG - still being considered for surgery' };

                                              
for iC = 1:numel(col_tgh)
    if ~isempty(col_tgh{iC})
        for iT = 1:size(col_tgh{iC},1)
            if isempty(dan_dta(strcmpi(dan_dta(:,iC),col_tgh{iC}{iT,2}),iC))
                dan_dta(strcmpi(dan_dta(:,iC),col_tgh{iC}{iT,2}),iC) = {col_tgh{iC}{strcmpi(col_tgh{iC}(:,1),'empty'),1}};
            else
                dan_dta(strcmpi(dan_dta(:,iC),col_tgh{iC}{iT,2}),iC) = {col_tgh{iC}{iT,1}};
            end
        end
    end
end

% [ dan_dta(:,strcmpi(dan_dta_col,'Surgery')) dan_dta(:,strcmpi(dan_dta_col,'SurgerySide'))]
% unique(dan_dta(:,strcmpi(dan_dta_col,'Surgery')))

%% Modify Rebecca Data
% [ reb_dta, reb_dta_sbj, reb_dta_col] = ejk_dta_frm( fcfg );

clear col_tgh

reb_dta     = reb_dta(:,[ strcmpi(reb_dta_col,'Type of Surgery - Comment') ]);
reb_dta_col = {'Surgery'};
reb_dta(:,end+1) = reb_dta(:,strcmpi(dan_dta_col,'Surgery'));
reb_dta_col      = [ reb_dta_col {'SurgerySide'} ];

col_tgh{find(strcmpi(reb_dta_col,'Surgery'))} = { 'ATL'  'Right ATL'  ; ...
                                                  'ATL'  'Right ATL (open) & repair of basal TL encephalocele'  ; ...
                                                  ...
                                                  'SLAH' 'Left SLAH' ; ...
                                                  'SLAH' 'Left SLAH (pt. developed PNES after surgery and has performed invalidly); had small left subdural bleed; seemed to be significantly impaired as well; likely had an unrecognized encephalitis; died of SUDEP during May of 2019' ; ...
                                                  'SLAH' 'Right SLAH' ; ...
                                                  'SLAH' 'Right SLAH (with DNET ablation attached to hippocampus)' ; ...
                                                  'SLAH' 'Right SLAH and ablation of a closely neighbording subcentimeter cavernous malformation' ; ...
                                                  ...
                                                  'Multiple' '1)  Left SLAH (plus - medial temporal pole ; parahippocampal gyrus; enthorhinal cortex ablation) 8/2018;  2) complete left ATL  02/2020' ; ...
                                                  'Multiple' '1) Left TL Gamma Knife Procedure (2007-UCSF) - failed to hit medial TL structures that were targeted. Insetead; caused necrotic change to the inferior longitudinal fasciculus; 2) Left SLAH (12/2016 - Emory); 3) SLA of left parahippocampal gyrus and residual amygdala (9/2017)' ; ...
                                                  'Multiple' '1) Patient had RF ablation of right amygdala; 2) Pt. had right SLAH' ; ...
                                                  'Multiple' 'Left SAH  (had left SLAH-but went directly to ATL without NP testing or neuroimaging in between)' ; ...
                                                  'Multiple' 'Left SLAH (10/2013) followed by further ablation with laser to clean up first on 1/2015)' ; ...
                                                  'Multiple' 'Left SLAH (followed by a second SLAH)' ; ...
                                                  'Multiple' 'Left SLAH (further clean up with SLAH)' ; ...
                                                  'Multiple' 'Left SLAH (with repeat SLAH for cleaning up residual tissue)' ; ...
                                                  'Multiple' 'Left TL cav. Mal. Ablation' ; ...
                                                  'Multiple' 'Prior right TL procedure; presented as bilteral TL onset now; and received bi-temporal RNS implantation' ; ...
                                                  'Multiple' 'Pt. had a prior left lateral TL resection and multiple subpial transection at MCG. Still had seizures and underwent SEEG at Emory' ; ...
                                                  'Multiple' 'Right SAH (open)  followed by right SAH procedure to clean up residual hippocampus' ; ...
                                                  'Multiple' 'Right SAH (open) (followed by right ATL procedure when not seizure free from first SAH)' ; ...
                                                  'Multiple' 'Right SLA (got parahippocampal gyrus only); followed by right SLAH; followed by Right ATL' ; ...
                                                  'Multiple' 'Right SLAH (later had second ablation to clean up first); Then had a right ATL' ; ...
                                                  'Multiple' 'Right amygdalotomy followed by limited right open resection after invasive monitoring' ; ...
                                                  'Multiple' 'limited right ATL resection; later followed with limited RF ablation following invasive monitoring' ; ...
                                                  'Multiple' 'prior right ATL; right inferior-posterior resection (extension)' ; ...
                                                  ...
                                                  'Other' 'Left SAH' ; ...
                                                  'Other' 'Left SAH (open)' ; ...
                                                  'Other' 'Left SLA' ; ...
                                                  'Other' 'Left SLA (involving the encephalocele in the left temporal pole)' ; ...
                                                  'Other' 'Left SLA - cavenous malformation (Location= temporal pole)' ; ...
                                                  'Other' 'Left SLA of DNET involving the posterior region of the middle temporal gyrus.' ; ...
                                                  'Other' 'Left SLA of anterior medial temporal pole cavernous malfromation' ; ...
                                                  'Other' 'Left SLA of cavernous malformation involving the white matter of the left fusiform gyrus' ; ...
                                                  'Other' 'Left SLA of left TL cavernous malformation (involving middle temporal gyrus)' ; ...
                                                  'Other' 'Left SLA of small DNET in medial TL region' ; ...
                                                  'Other' 'Left inferior-posterior TL resection from grid (visual word-form area)' ; ...
                                                  'Other' 'Right Frontal lobe open resection (involved right frontal operculum)' ; ...
                                                  'Other' 'Right SAH' ; ...
                                                  'Other' 'Right SAH (open)' ; ...
                                                  'Other' 'Right SAH (open);' ; ...
                                                  'Other' 'Right SLA of cavernous malformation in the right anterior temporal pole' ; ...
                                                  'Other' 'Right TL cav mal ablation' ; ...
                                                  'Other' 'Right TL open resection (not standard)' ; ...
                                                  'Other' 'Right insulotomy (open)' ; ...
                                                  'Other' 'SEEG/Followed by SLA of right temporal pole extending to amygdala but sparing hippocampus' ; ...
                                                  'Other' 'SLA of cavernous malformation involving the left posterior inferior temporal gyrus' ; ...
                                                  'Other' 'SLA of cavernous malformation of the right anterior temporal pole' ; ...
                                                  'Other' 'SLA of left TL tumor (probable DNET)' ; ...
                                                  'Other' 'SLA of left temporal polar and pterional cortex and left amygdala in a patient with an encephalocele associated with epilepsy' ; ...
                                                  'Other' 'SLA of right anterior TL (more than a SLAH/almost like an ATL done with laser)' ; ...
                                                   ...
                                                  'None' 'empty' };


                                              
for iC = 1:numel(col_tgh)
    if ~isempty(col_tgh{iC})
        for iT = 1:size(col_tgh{iC},1)
            if strcmpi(col_tgh{iC}{iT,2},'empty')
                reb_dta(cellfun(@isempty,reb_dta(:,strcmpi(reb_dta_col,'Surgery'))),iC) = {col_tgh{iC}{iT,1}};
            else
                reb_dta(strcmpi(reb_dta(:,iC),col_tgh{iC}{iT,2}),iC) = {col_tgh{iC}{iT,1}};
            end
        end
    end
end

% [ reb_dta(:,strcmpi(reb_dta_col,'Surgery')) reb_dta(:,strcmpi(reb_dta_col,'SurgerySide'))]
% unique(reb_dta(:,strcmpi(reb_dta_col,'Surgery')))                                              
                                              
%% Combine Emory
tot_sbj = unique( [dan_dta_sbj ; reb_dta_sbj] );

out_emy_dta = cell(numel(tot_sbj),5);
for iS = 1:numel(tot_sbj)
    
    out_emy_dta{iS,1} = tot_sbj{iS};
    
    dan_ind = find(strcmpi(dan_dta_sbj,tot_sbj{iS}));
    reb_ind = find(strcmpi(reb_dta_sbj,tot_sbj{iS}));
    
    if ~isempty(dan_ind)
        out_emy_dta{iS,2} = dan_dta{dan_ind,strcmpi(dan_dta_col,'Surgery')};
        out_emy_dta{iS,4} = dan_dta{dan_ind,strcmpi(dan_dta_col,'SurgerySide')};
    else
        out_emy_dta{iS,2} = '';
        out_emy_dta{iS,4} = '';
    end
    
    if ~isempty(reb_ind)
        out_emy_dta{iS,3} = reb_dta{reb_ind,strcmpi(reb_dta_col,'Surgery')};
        out_emy_dta{iS,5} = reb_dta{reb_ind,strcmpi(reb_dta_col,'SurgerySide')};
    else
        out_emy_dta{iS,3} = '';
        out_emy_dta{iS,5} = '';
    end
    
end

cell2csv([ out_dir '/' 'Emory_Check.csv'],[ {'sbj_nme' 'DanSurgery' 'RebeccaSurgery' 'DanNotes' 'RebeccaNotes' } ; out_emy_dta ])

%% Load UCSD
fcfg = [];
fcfg.red_fle = [ dta_dir '/' 'Redcap_2022_07_22.csv'];
fcfg.sep     = '|';
[~ , ~ , ~ , fst_sbj_cog, ~, fst_sbj_srg] = mmil_load_redcap(fcfg);

fcfg = [];
fcfg.red_fle = [ dta_dir '/' 'Redcap_2023_01_12.csv'];
fcfg.sep     = '|';
[~ , ~ , ~ , scd_sbj_cog, ~, scd_sbj_srg] = mmil_load_redcap(fcfg);

%% Combine UCSD
san_dgo_out = cell(numel(scd_sbj_cog.sbj_nme),7);

for iS = 1:numel(scd_sbj_cog.sbj_nme)
    
    san_dgo_out{iS,1} = scd_sbj_cog.sbj_nme{iS};
    
    fst_ind = find(strcmpi(fst_sbj_cog.sbj_nme,scd_sbj_cog.sbj_nme{iS}));
    
    if ~isempty(fst_ind)
        san_dgo_out{iS,2} = fst_sbj_srg.srg_typ{fst_ind};
        san_dgo_out{iS,3} = fst_sbj_cog.log_mem_raw_scr_two(fst_ind);
        san_dgo_out{iS,4} = fst_sbj_cog.log_mem_raw_scr_two_pst(fst_ind);
    else
        san_dgo_out{iS,2} = 'MISSING';
        san_dgo_out{iS,3} = NaN;
        san_dgo_out{iS,4} = NaN;
    end
    
    san_dgo_out{iS,5} = scd_sbj_srg.srg_typ{iS};
    san_dgo_out{iS,6} = scd_sbj_cog.log_mem_raw_scr_two(iS);
    san_dgo_out{iS,7} = scd_sbj_cog.log_mem_raw_scr_two_pst(iS);
    
end

cell2csv([ out_dir '/' 'UCSD_Check.csv'],[ {'sbj_nme' '1st_Surgery' '1st_LM_pre' '1st_LM_post' '2nd_Surgery' '2nd_LM_pre' '2nd_LM_post'} ; san_dgo_out ])

%% Combine w/ existing data
fcfg = [];
fcfg.dta_loc = [ out_dir '/' 'FinalSample.csv'];
fcfg.dta_col = 2;
[ fnl_dta, fnl_dta_sbj, fnl_dta_col] = ejk_dta_frm( fcfg );

load([ out_dir '/' 'grp_img_qal.mat'])
grp_inc = fnl_dta_sbj( [ grp.surgery.pst_cog_dti.ltle_slah ; grp.surgery.pst_cog_dti.ltle_atl ; grp.surgery.pst_cog_dti.rtle_slah ; grp.surgery.pst_cog_dti.rtle_atl ]);

fcfg = [];
fcfg.dta_loc = [ out_dir '/' 'Emory_Check.csv'];
fcfg.dta_col = 2;
[ emy_dta, emy_dta_sbj, emy_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ out_dir '/' 'UCSD_Check.csv'];
fcfg.dta_col = 2;
[ sdg_dta, sdg_dta_sbj, sdg_dta_col] = ejk_dta_frm( fcfg );

tot_sbj = unique( [ fnl_dta_sbj ; emy_dta_sbj ; sdg_dta_sbj] );

out_tbl = cell(numel(tot_sbj),9);
for iS = 1:numel(tot_sbj)
    
    fnl_ind = strcmpi(fnl_dta_sbj,tot_sbj{iS});
    emy_ind = strcmpi(emy_dta_sbj,tot_sbj{iS});
    sdg_ind = strcmpi(sdg_dta_sbj,tot_sbj{iS});
    grp_ind = strcmpi(grp_inc,tot_sbj{iS});
    
    out_tbl{iS,1} = tot_sbj{iS};
    
    if sum(grp_ind)==0
        out_tbl{iS,2} = '';
    else
        out_tbl{iS,2} = 'INCLUDED';
    end
    
    if sum(fnl_ind)==0
        out_tbl{iS,3} = '';
        out_tbl{iS,4} = [];
    else
        out_tbl{iS,3} = fnl_dta{fnl_ind,strcmpi(fnl_dta_col,'SurgeryType')};
        out_tbl{iS,4} = fnl_dta{fnl_ind,strcmpi(fnl_dta_col,'lm2_chg')};
    end
    
    if sum(sdg_ind)==0
        out_tbl{iS,5} = '';
        out_tbl{iS,6} = [];
        out_tbl{iS,7} = [];
    else
        out_tbl{iS,5} = sdg_dta{sdg_ind,strcmpi(sdg_dta_col,'1st_Surgery')};
        out_tbl{iS,6} = sdg_dta{sdg_ind,strcmpi(sdg_dta_col,'1st_LM_post')};
        out_tbl{iS,7} = sdg_dta{sdg_ind,strcmpi(sdg_dta_col,'2nd_LM_post')};
    end
    
    if sum(emy_ind)==0
        out_tbl{iS,8} = '';
        out_tbl{iS,9} = '';
    else
        out_tbl{iS,8} = emy_dta{emy_ind,strcmpi(emy_dta_col,'DanSurgery')};
        out_tbl{iS,9} = emy_dta{emy_ind,strcmpi(emy_dta_col,'RebeccaSurgery')};
    end
    
end

cell2csv([ out_dir '/' 'Total_Check.csv'],[ {'sbj_nme' 'Included' 'CurrentSurgery' 'CurrentLMpost' 'RedcapSurgery' 'Redcap1stLMpost' 'Redcap2ndLMpost' 'DanSurgery' 'RebeccaSurgery'} ; out_tbl ])














%% Load data
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

%% Define divisions of interest
usd_ind = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'Location')),'UCSD'));
usf_ind = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'Location')),'UCSF'));
emy_ind = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'Location')),'Emory'));

has_pre_cog = find(~isnan(cell2mat(cog_dta(:,strcmpi(cog_dta_col,'lm2_pre')))) | ...
                   ~isnan(cell2mat(cog_dta(:,strcmpi(cog_dta_col,'vp2_pre')))) | ...
                   ~isnan(cell2mat(cog_dta(:,strcmpi(cog_dta_col,'cv2_pre')))) );
has_pst_cog = find(~isnan(cell2mat(cog_dta(:,strcmpi(cog_dta_col,'lm2_pst')))) | ...
                   ~isnan(cell2mat(cog_dta(:,strcmpi(cog_dta_col,'vp2_pst')))) | ...
                   ~isnan(cell2mat(cog_dta(:,strcmpi(cog_dta_col,'cv2_pst')))) ); 
               
ltl_tle = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'SideOfSeizureFocus')),'L'));
rgh_tle = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'SideOfSeizureFocus')),'R'));
btl_tle = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'SideOfSeizureFocus')),'B'));
con_trl     = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'Location')),'HC'));

slh_tle = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'SurgeryType')),'SLAH'));
atl_tle = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'SurgeryType')),'ATL'));

lm2_wms_thr = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'wms_lm2')),'III'));
lm2_wms_for = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'wms_lm2')),'IV'));

vp2_wms_two = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'wms_vp2')),'II'));
vp2_wms_thr = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'wms_vp2')),'III'));
vp2_wms_for = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'wms_vp2')),'IV'));

%% Grouping 1: Site
grp.site.control = con_trl;
grp.site.emory   = emy_ind;
grp.site.ucsf    = usf_ind;
grp.site.ucsd    = usd_ind;

%% Grouping 1.1: Diagnosis 
grp.diagnosis.control = con_trl;
grp.diagnosis.ltle   = ltl_tle;
grp.diagnosis.rtle   = rgh_tle;

%% Grouping 2: presurgical-x-condition
grp.pre_cog.control = intersect(has_pre_cog,con_trl);
grp.pre_cog.ltle   = intersect(has_pre_cog,ltl_tle);
grp.pre_cog.rtle    = intersect(has_pre_cog,rgh_tle);

%% Grouping 2.1: Total Sample Diagnosis 
grp.total_pre_cog.tle_hc = [ grp.pre_cog.control ; grp.pre_cog.ltle ; grp.pre_cog.rtle];

%% Grouping 3: postsurgical-x-condition
grp.pst_cog.ltle    = intersect(has_pst_cog,ltl_tle);
grp.pst_cog.rtle    = intersect(has_pst_cog,rgh_tle);

%% Grouping 4: presurgical-x-condition-x-imaging
% use_mri_dta     = find( ~isnan(cell2mat(mri_dev_dta(:,strcmpi(mri_dev_dta_col,'Left_Hippocampus')))) );
% use_dti_fib_dta = find( ~isnan(cell2mat(mri_dev_dta(:,strcmpi(dti_fib_dev_dta_col,'L_IFO')))) );
% use_dti_wmp_dta = find( ~isnan(cell2mat(mri_dev_dta(:,strcmpi(dti_wmp_dev_dta_col,'lh_parsorbitalis')))) );
% 
% grp.pre_cog_img.control = intersect(intersect(has_pre_cog,con_trl),use_dti_fib_dta);
% grp.pre_cog_img.ltle    = intersect(intersect(has_pre_cog,ltl_tle),use_dti_fib_dta);
% grp.pre_cog_img.rtle    = intersect(intersect(has_pre_cog,rgh_tle),use_dti_fib_dta);

%% Grouping 5: postsurgical-x-condition-x-imaging
% use_mri_dta     = find( ~isnan(cell2mat(mri_dev_dta(:,strcmpi(mri_dev_dta_col,'Left_Hippocampus')))) );
% use_dti_fib_dta = find( ~isnan(cell2mat(mri_dev_dta(:,strcmpi(dti_fib_dev_dta_col,'L_IFO')))) );
% use_dti_wmp_dta = find( ~isnan(cell2mat(mri_dev_dta(:,strcmpi(dti_wmp_dev_dta_col,'lh_parsorbitalis')))) ); 
% 
% grp.pst_cog_img.ltle    = intersect(intersect(has_pst_cog,ltl_tle),use_dti_fib_dta);
% grp.pst_cog_img.rtle    = intersect(intersect(has_pst_cog,rgh_tle),use_dti_fib_dta);

%% Grouping 6: SurgeryType
grp.surgery.ltle_slah = intersect(ltl_tle,slh_tle);
grp.surgery.ltle_atl  = intersect(ltl_tle,atl_tle);
grp.surgery.rtle_slah = intersect(rgh_tle,slh_tle);
grp.surgery.rtle_atl  = intersect(rgh_tle,atl_tle);

%% Grouping 7: WMS version - explore
% grp.wms_ver.lm2_two   = intersect(find(strcmpi(wms_hld(:,2),'II')),grp.site.emory);
% grp.wms_ver.lm2_three = intersect(find(strcmpi(wms_hld(:,2),'III')),grp.site.emory);
% grp.wms_ver.lm2_four  = intersect(find(strcmpi(wms_hld(:,2),'IV')),grp.site.emory);
% grp.wms_ver.vp2_two   = intersect(find(strcmpi(wms_hld(:,3),'II')),grp.site.emory);
% grp.wms_ver.vp2_three = intersect(find(strcmpi(wms_hld(:,3),'III')),grp.site.emory);
% grp.wms_ver.vp2_four  = intersect(find(strcmpi(wms_hld(:,3),'IV')),grp.site.emory);

%% Grouping 8: L-TLE only
grp.surgery_ltle.ltle_slah = intersect(ltl_tle,slh_tle);
grp.surgery_ltle.ltle_atl  = intersect(ltl_tle,atl_tle);

%% Grouping 9: Surgery-Only
grp.surgery_only.slah = slh_tle;
grp.surgery_only.atl  = atl_tle;

%% Grouping 10: LM2 WMS version
grp.lm2_wms_ver.III = lm2_wms_thr;
grp.lm2_wms_ver.IV  = lm2_wms_for;

%% Grouping 11: VP2 WMS version 
grp.vp2_wms_ver.II  = vp2_wms_two;
grp.vp2_wms_ver.III = vp2_wms_thr;
grp.vp2_wms_ver.IV  = vp2_wms_for;

%% Make group names
grp_nme = fieldnames(grp);
for iG = 1:numel(grp_nme)
    grp_nme_sub = fieldnames(grp.(grp_nme{iG}));
    grp.grp_cat.(grp_nme{iG}) = repmat({''},numel(cln_dta_sbj),1);
    for iGS = 1:numel(grp_nme_sub)
        grp.grp_cat.(grp_nme{iG})(grp.(grp_nme{iG}).(grp_nme_sub{iGS})) = {grp_nme_sub{iGS}};
    end
end

%% Save out
save([ prj_dir '/' prj_nme '/' 'Data' '/' 'grp.mat'],'grp');

%% Display groups
grp_nme = fieldnames(grp);
dsp_txt = '';
for iG = 1:numel(grp_nme)-1
    dsp_txt = [dsp_txt '\n\n' ' ' grp_nme{iG}];
    grp_nme_sub = fieldnames(grp.(grp_nme{iG}));
    for iGS = 1:numel(grp_nme_sub)
        dsp_txt = [dsp_txt '\n' '   ' grp_nme_sub{iGS} ':' ' ' num2str(numel(grp.(grp_nme{iG}).(grp_nme_sub{iGS})))];
    end
end
fprintf([dsp_txt '\n\n\n'])








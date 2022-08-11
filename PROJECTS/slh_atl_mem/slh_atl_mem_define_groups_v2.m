%% Load data
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'fiber_FA_dev.csv'];
fcfg.dta_col = 2;
[ dti_dta, dti_dta_sbj, dti_dta_col] = ejk_dta_frm( fcfg );

%% Define divisions of interest
% Site indices
usd_ind = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'Location')),'UCSD'));
usf_ind = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'Location')),'UCSF'));
emy_ind = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'Location')),'Emory'));

% Data indices
has_pre_cog = find(~isnan(cell2mat(cog_dta(:,strcmpi(cog_dta_col,'lm2_pre')))) | ...
                   ~isnan(cell2mat(cog_dta(:,strcmpi(cog_dta_col,'vp2_pre')))) | ...
                   ~isnan(cell2mat(cog_dta(:,strcmpi(cog_dta_col,'cv2_pre')))) );
has_pst_cog = find(~isnan(cell2mat(cog_dta(:,strcmpi(cog_dta_col,'lm2_pst')))) | ...
                   ~isnan(cell2mat(cog_dta(:,strcmpi(cog_dta_col,'vp2_pst')))) | ...
                   ~isnan(cell2mat(cog_dta(:,strcmpi(cog_dta_col,'cv2_pst')))) );                
has_pre_dti = find(~isnan(cell2mat(dti_dta(:,strcmpi(dti_dta_col,'L_Unc')))) | ...
                   ~isnan(cell2mat(dti_dta(:,strcmpi(dti_dta_col,'R_IFO')))) | ...
                   ~isnan(cell2mat(dti_dta(:,strcmpi(dti_dta_col,'L_Fx')))) );

% Disease indices
ltl_tle = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'SideOfSeizureFocus')),'L'));
rgh_tle = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'SideOfSeizureFocus')),'R'));
btl_tle = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'SideOfSeizureFocus')),'B'));
con_trl     = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'Location')),'HC'));

% Surgery indices
slh_tle = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'SurgeryType')),'SLAH'));
atl_tle = find(strcmpi(cln_dta(:,strcmpi(cln_dta_col,'SurgeryType')),'ATL'));

%% Grouping 1: Site
grp.site.pre_cog.control = intersect(con_trl,has_pre_cog);
grp.site.pre_cog.emory   = intersect(emy_ind,has_pre_cog);
grp.site.pre_cog.ucsf    = intersect(usf_ind,has_pre_cog);
grp.site.pre_cog.ucsd    = intersect(usd_ind,has_pre_cog);

grp.site.pre_cog_dti.control = intersect(intersect(con_trl,has_pre_cog),has_pre_dti);
grp.site.pre_cog_dti.emory   = intersect(intersect(emy_ind,has_pre_cog),has_pre_dti);
grp.site.pre_cog_dti.ucsf    = intersect(intersect(usf_ind,has_pre_cog),has_pre_dti);
grp.site.pre_cog_dti.ucsd    = intersect(intersect(usd_ind,has_pre_cog),has_pre_dti);

grp.site.pst_cog.control = intersect(con_trl,has_pst_cog);
grp.site.pst_cog.emory   = intersect(emy_ind,has_pst_cog);
grp.site.pst_cog.ucsf    = intersect(usf_ind,has_pst_cog);
grp.site.pst_cog.ucsd    = intersect(usd_ind,has_pst_cog);

grp.site.pst_cog_dti.control = intersect(intersect(con_trl,has_pst_cog),has_pre_dti);
grp.site.pst_cog_dti.emory   = intersect(intersect(emy_ind,has_pst_cog),has_pre_dti);
grp.site.pst_cog_dti.ucsf    = intersect(intersect(usf_ind,has_pst_cog),has_pre_dti);
grp.site.pst_cog_dti.ucsd    = intersect(intersect(usd_ind,has_pst_cog),has_pre_dti);

%% Grouping 2: Pre-Surgical by Diagnosis 
grp.diagnosis.pre_cog.control = intersect(con_trl,has_pre_cog);
grp.diagnosis.pre_cog.ltle    = intersect(ltl_tle,has_pre_cog);
grp.diagnosis.pre_cog.rtle    = intersect(rgh_tle,has_pre_cog);

grp.diagnosis.pre_cog_dti.control = intersect(intersect(con_trl,has_pre_cog),has_pre_dti);
grp.diagnosis.pre_cog_dti.ltle    = intersect(intersect(ltl_tle,has_pre_cog),has_pre_dti);
grp.diagnosis.pre_cog_dti.rtle    = intersect(intersect(rgh_tle,has_pre_cog),has_pre_dti);

%% Grouping 3: Total Sample Diagnosis 
grp.total_pre_cog.pre_cog.tle_hc     = [ con_trl ; ltl_tle ; rgh_tle ];
grp.total_pre_cog.pre_cog_dti.tle_hc = intersect(grp.total_pre_cog.pre_cog.tle_hc, has_pre_dti);

%% Grouping 4: Post-Surgical by Surgery 
grp.surgery.pst_cog.ltle_slah = intersect(intersect(ltl_tle,slh_tle),has_pst_cog);
grp.surgery.pst_cog.ltle_atl  = intersect(intersect(ltl_tle,atl_tle),has_pst_cog);
grp.surgery.pst_cog.rtle_slah = intersect(intersect(rgh_tle,slh_tle),has_pst_cog);
grp.surgery.pst_cog.rtle_atl  = intersect(intersect(rgh_tle,atl_tle),has_pst_cog);

grp.surgery.pst_cog_dti.ltle_slah = intersect(intersect(intersect(ltl_tle,slh_tle),has_pst_cog), has_pre_dti);
grp.surgery.pst_cog_dti.ltle_atl  = intersect(intersect(intersect(ltl_tle,atl_tle),has_pst_cog), has_pre_dti);
grp.surgery.pst_cog_dti.rtle_slah = intersect(intersect(intersect(rgh_tle,slh_tle),has_pst_cog), has_pre_dti);
grp.surgery.pst_cog_dti.rtle_atl  = intersect(intersect(intersect(rgh_tle,atl_tle),has_pst_cog), has_pre_dti);

%% Make group names
% grp_nme = fieldnames(grp);
% for iG = 1:numel(grp_nme)
%     grp_nme_sub = fieldnames(grp.(grp_nme{iG}));
%     for iGS = 1:numel(grp_nme_sub)
%         cat_nme = fieldnames(grp.(grp_nme{iG}).(grp_nme_sub{iGS}));
%         grp.grp_cat.(grp_nme{iG}).(grp_nme_sub{iGS}) = repmat({''},numel(cln_dta_sbj),1);
%         for iCA = 1:numel(cat_nme)
%             grp.grp_cat.(grp_nme{iG}).(grp_nme_sub{iGS})(grp.(grp_nme{iG}).(grp_nme_sub{iGS}).(cat_nme{iCA})) = {cat_nme{iCA}};
%         end
%     end
% end

%% Save out
save([ prj_dir '/' prj_nme '/' 'Data' '/' 'grp.mat'],'grp');

%% Display groups
% grp_nme = fieldnames(grp);
% dsp_txt = '';
% for iG = 1:numel(grp_nme)-1
%     dsp_txt = [dsp_txt '\n\n' ' ' grp_nme{iG}];
%     grp_nme_sub = fieldnames(grp.(grp_nme{iG}));
%     for iGS = 1:numel(grp_nme_sub)
%         dsp_txt = [dsp_txt '\n' '   ' grp_nme_sub{iGS} ':' ' ' num2str(numel(grp.(grp_nme{iG}).(grp_nme_sub{iGS})))];
%     end
% end
% fprintf([dsp_txt '\n\n\n'])








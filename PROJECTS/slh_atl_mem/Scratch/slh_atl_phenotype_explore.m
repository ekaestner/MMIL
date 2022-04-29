load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )
grp_nme = fieldnames(grp);

dta_fld = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/Data/';

out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/slh_atl_mem/explore/phenotyping';

cln_dta = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv' ]);
cln_dta_col = cln_dta(1,2:end);
cln_dta_sbj = cln_dta(2:end,1);
cln_dta     = cln_dta(2:end,2:end);

cog_dta     = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv']);
cog_dta_col = cog_dta(1,2:end);
cog_dta_sbj = cog_dta(2:end,1);
cog_dta     = cog_dta(2:end,2:end);

%% %%%%%%%%%%%%%%%%%%%%%
% dta_sub_mse = { 'subcort_vol_238_LateralityIndex'              'subcort_vol_dev_LateralityIndex' ...
%                 'wmparc_FA_wm_aparc_annot_238_LateralityIndex' 'wmparc_FA_wm_aparc_annot_dev_LateralityIndex' ...
%                 'fiber_FA_238_LateralityIndex' 'fiber_FA_dev_LateralityIndex' };
% 
% roi_sub_fle = { ''               'fiber_FA_238_LateralityIndex_QC.csv' '' 'wmparc_FA_wm_aparc_annot_238_LateralityIndex_QC.csv' };
% roi_sub_mse = { 'cort_thick_ctx' 'fiber_FA'            'subcort_vol' 'wmparc_FA_wm' };
% roi_sub_nme = { { 'entorhinal' 'parahippocampal' } ...
%                 { 'Unc' 'ILF' } ...
%                 { 'hippocampus' 'Thalamus_Proper' } ...
%                 { 'entorhinal' 'parahippocampal' } };

%% Scatter
iG = 8; % 8 / 

cog_dta_one = { 'vp2_nor_scr_pst' 'log_mem_nor_scr_two_pst' 'vp2_nor_scr_pst' 'log_mem_nor_scr_two_pst' }; % 
cog_dta_two = { 'vp2_nor_scr_zsc' 'log_mem_nor_scr_two_zsc' 'vp2_nor_scr'     'log_mem_nor_scr_two' }; % 
iC = 1;

fcfg = [];

fcfg.xdt     = { cell2mat(cog_dta( grp.(grp_nme{iG}) ,strcmpi(cog_dta_col,cog_dta_one{iC}) )) };
fcfg.ydt     = { cell2mat(cog_dta( grp.(grp_nme{iG}) ,strcmpi(cog_dta_col,cog_dta_two{iC}) )) };

fcfg.fce_col = { rgb('deep green') };
fcfg.edg_col = { rgb('black') };

fcfg.xlb = { cog_dta_one{iC} };
fcfg.ylb = { cog_dta_two{iC}  };

fcfg.hln     = -1.28;
fcfg.hln_col = rgb('red');

fcfg.vln     = -1.28;
fcfg.vln_col = rgb('red');

fcfg.xlm = [ -5 5 ];
fcfg.ylm = [ -5 1 ];

fcfg.out_dir = out_dir;
fcfg.out_nme = [ grp_nme{iG} '_' cog_dta_one{iC} '_BY_' cog_dta_two{iC} ];

ejk_scatter(fcfg)
            
%% Table
cog_dta_one = { 'log_mem_nor_scr_two_cat'     'vp2_nor_scr_cat' }; % 
cog_dta_two = { 'log_mem_nor_scr_two_pst_cat' 'vp2_nor_scr_pst_cat' }; % 
iC = 2;

grp_nme = fieldnames(grp);
iG = 8; % 8 / 

pre_hld = cog_dta( grp.(grp_nme{iG}) ,strcmpi(cog_dta_col,cog_dta_one{iC}) );
    pre_hld(cellfun(@(x) x==1,cellfun(@isnumeric,pre_hld,'uni',0))) = {''};
pst_hld = cog_dta( grp.(grp_nme{iG}) ,strcmpi(cog_dta_col,cog_dta_two{iC}) );
    pst_hld(cellfun(@(x) x==1,cellfun(@isnumeric,pst_hld,'uni',0))) = {''};

[ crs_tbl, ~, ~, crs_lbl ] = crosstab( categorical(pre_hld), categorical(pst_hld));

[ {''} strcat('POST_',crs_lbl(:,2))' ; strcat('PRE_',crs_lbl(:,1)) num2cell(crs_tbl) ] 



[ cog_dta( grp.(grp_nme{iG}), 8) cog_dta( grp.(grp_nme{iG}), 10) ]
[ cog_dta( grp.(grp_nme{iG}), 7) cog_dta( grp.(grp_nme{iG}), 9) ]

%% Color Scatter
grp_nme = fieldnames(grp);
iG = 8; % 8 / 

cog_dta_one = { 'vp2_nor_scr_pst' 'log_mem_nor_scr_two_pst' }; % 
cog_dta_two = { 'vp2_nor_scr_zsc' 'log_mem_nor_scr_two_zsc' }; % 
iC = 1;

iR = 4;
iRO = 1;
fcfg = [];
fcfg.dta_loc = [ dta_fld '/' roi_sub_fle{iR}];
fcfg.dta_col = 5;
[ neu_dta, neu_dta_sbj, neu_dta_col] = ejk_dta_frm( fcfg );

top_pct = 1;
cfg.col_map = { rgb('bright red') rgb('dark purple') rgb('bright blue') };
col_map = [];
for iC = 1:numel(cfg.col_map)-1
    col_map = [col_map ; [linspace(cfg.col_map{iC}(1),cfg.col_map{iC+1}(1),ceil(1000*top_pct/(numel(cfg.col_map)-1)))' linspace(cfg.col_map{iC}(2),cfg.col_map{iC+1}(2),ceil(1000*top_pct/(numel(cfg.col_map)-1)))' linspace(cfg.col_map{iC}(3),cfg.col_map{iC+1}(3),ceil(1000*top_pct/(numel(cfg.col_map)-1)))']; ];
end

neu_dta_use = cell2mat(neu_dta(grp.(grp_nme{iG}),strcmpi(neu_dta_col,roi_sub_nme{iR}{iRO})));
xdt_dta_use = cell2mat(cog_dta( grp.(grp_nme{iG}) ,strcmpi(cog_dta_col,cog_dta_one{iC}) ));
ydt_dta_use = cell2mat(cog_dta( grp.(grp_nme{iG}) ,strcmpi(cog_dta_col,cog_dta_two{iC}) ));
min_val = min(neu_dta_use);
max_val = max(neu_dta_use);
min_val = min_val - ((max_val-min_val)*.05);   

clear xdt ydt fce_col
for iCL = 1:numel(neu_dta_use)    
    xdt{iCL} = xdt_dta_use(iCL);
    ydt{iCL} = ydt_dta_use(iCL);
    if ~isnan(neu_dta_use(iCL))
    fce_col{iCL} = col_map(ceil(((neu_dta_use(iCL)-min_val) / (max_val - min_val))*1000)-1,:);
    else
        fce_col{iCL} = rgb('white');
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.xdt     = xdt;
fcfg.ydt     = ydt;

fcfg.fce_col = fce_col;
fcfg.edg_col = repmat({rgb('black')},1,numel(xdt));

fcfg.xlb = { cog_dta_one{iC} };
fcfg.ylb = { cog_dta_two{iC}  };

fcfg.hln     = -1.28;
fcfg.hln_col = rgb('red');

fcfg.vln     = -1.28;
fcfg.vln_col = rgb('red');

fcfg.jtr = 0;
fcfg.xlm = [ -5 5 ];
fcfg.ylm = [ -5 1 ];

fcfg.out_dir = out_dir;
fcfg.out_nme = [ grp_nme{iG} '_' cog_dta_one{iC} '_BY_' cog_dta_two{iC} '_' roi_sub_nme{iR}{iRO} ];

ejk_scatter(fcfg)

%% Whole Group
ejk_chk_dir([ out_dir '/' 'post_scatters' '/' ])
ejk_chk_dir([ out_dir '/' 'pre_scatters' '/' ])
ejk_chk_dir([ out_dir '/' 'phenotype_scatters_238' '/' ])
ejk_chk_dir([ out_dir '/' 'phenotype_scatters_dev' '/' ])

cog_dta_one = { 'log_mem_nor_scr_two_pst' }; % 
cog_dta_two = { 'log_mem_nor_scr_two_zsc' }; % 

roi_sub_fle = { 'subcort_vol_XREP_norm_IntracranialVolume_QC.csv' ...
                'subcort_vol_XREP_LateralityIndex_QC.csv' ...
                'fiber_FA_XREP_QC.csv' ...
                'fiber_FA_XREP_LateralityIndex_QC.csv' ...
                'wmparc_FA_wm_aparc_annot_XREP_QC.csv' ...
                'wmparc_FA_wm_aparc_annot_XREP_LateralityIndex_QC.csv' };
roi_sub_mse = { 'volume'       ...
                'volume_LI'    ...
                'fiber_FA'     ...
                'fiber_FA_LI'  ...
                'wmparc_FA_wm' ...
                'wmparc_FA_wm_LI' };
roi_sub_nme = { { 'Left_Hippocampus' } ...
                { 'Hippocampus' } ...
                { 'L_Unc' 'L_ILF' } ...
                { 'Unc' 'ILF' } ...                
                { 'lh_entorhinal' 'lh_lateralorbitofrontal' 'lh_medialorbitofrontal' 'lh_parsopercularis' } ...
                { 'entorhinal' 'lateralorbitofrontal' 'medialorbitofrontal' 'parsopercularis'  } };
  
fcfg = [];
fcfg.dta_loc = [ dta_fld '/' 'Cognitive_QC.csv' ];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

for iD = 1:numel(roi_sub_fle)
    
    % Load
    fcfg = [];
    fcfg.dta_loc = [ dta_fld '/' strrep(roi_sub_fle{iD},'XREP','238') ];
    fcfg.dta_col = 5;
    [ neu_238_dta, neu_238_dta_sbj, neu_238_dta_col] = ejk_dta_frm( fcfg );

    fcfg = [];
    fcfg.dta_loc = [dta_fld '/' strrep(roi_sub_fle{iD},'XREP','dev') ];
    fcfg.dta_col = 5;
    [ neu_dev_dta, neu_dev_dta_sbj, neu_dev_dta_col] = ejk_dta_frm( fcfg );    
    
    for iR = 1:numel(roi_sub_nme{iD})
        
        neu_238_col = find(strcmpi(neu_238_dta_col,roi_sub_nme{iD}{iR}));
        neu_dev_col = find(strcmpi(neu_dev_dta_col,roi_sub_nme{iD}{iR}));
               
        % Pre-operative (all)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cog_col_two = find(strcmpi(cog_dta_col,cog_dta_two{1}));
        
        fcfg = [];
        
        fcfg.xdt     = { cell2mat(neu_238_dta( grp.tle_controls_pre_3T_allSurg_all , neu_238_col )) ...
                         cell2mat(neu_dev_dta( grp.tle_controls_pre_3T_allSurg_all , neu_dev_col )) };
        fcfg.ydt     = { cell2mat(cog_dta( grp.tle_controls_pre_3T_allSurg_all, cog_col_two )) ...
                         cell2mat(cog_dta( grp.tle_controls_pre_3T_allSurg_all, cog_col_two )) };
        
        fcfg.fce_col = { rgb('deep green') rgb('burnt orange') };
        fcfg.edg_col = { rgb('black')      rgb('black') };
        
        fcfg.trd_lne = [ 1 1 ];
        
        fcfg.ylb = { mmil_spec_char(cog_dta_two{1},{'_'},{' '}) };
        fcfg.xlb = { mmil_spec_char(neu_238_dta_col{neu_238_col},{'_'},{' '})  };
        
        fcfg.hln     = -1.28;
        fcfg.hln_col = rgb('red');
        
        fcfg.vln     = -1.28;
        fcfg.vln_col = rgb('red');
        
        fcfg.ylm = [ -5 1 ];
        
        fcfg.out_dir = [ out_dir '/' 'pre_scatters' '/' ];
        fcfg.out_nme = [ roi_sub_mse{iD} '_' neu_238_dta_col{neu_238_col} '_BY_' cog_dta_two{1} ];
        
        ejk_scatter(fcfg)        
        
        % Post-operative scatter (LTLE)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cog_col_one = find(strcmpi(cog_dta_col,cog_dta_one{1}));
        
        fcfg = [];
        
        fcfg.xdt     = { cell2mat(neu_238_dta( grp.tle_post_3T_ATLonly_left , neu_238_col )) ...
                         cell2mat(neu_dev_dta( grp.tle_post_3T_ATLonly_left , neu_dev_col )) };
        fcfg.ydt     = { cell2mat(cog_dta( grp.tle_post_3T_ATLonly_left, cog_col_one )) ...
                         cell2mat(cog_dta( grp.tle_post_3T_ATLonly_left, cog_col_one )) };
        
        fcfg.fce_col = { rgb('deep green') rgb('burnt orange') };
        fcfg.edg_col = { rgb('black')      rgb('black') };
        
        fcfg.trd_lne = [ 1 1 ];
        
        fcfg.ylb = { mmil_spec_char(cog_dta_one{1},{'_'},{' '}) };
        fcfg.xlb = { mmil_spec_char(neu_238_dta_col{neu_238_col},{'_'},{' '})  };
        
        fcfg.hln     = -1.28;
        fcfg.hln_col = rgb('red');
        
        fcfg.vln     = -1.28;
        fcfg.vln_col = rgb('red');
        
        fcfg.ylm = [ -5 1 ];
        
        fcfg.out_dir = [ out_dir '/' 'post_scatters' '/' ];
        fcfg.out_nme = [ roi_sub_mse{iD} '_' neu_238_dta_col{neu_238_col} '_BY_' cog_dta_one{1} ];
        
        ejk_scatter(fcfg)        
        
        % Post-operative Phenotype Scatter (LTLE) - DEV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make color map
        cfg.col_map = { rgb('neon red') rgb('red') rgb('dark purple') rgb('blue') rgb('neon blue') };
        col_map = [];
        for iC = 1:numel(cfg.col_map)-1
            col_map = [col_map ; [linspace(cfg.col_map{iC}(1),cfg.col_map{iC+1}(1),ceil(1000/(numel(cfg.col_map)-1)))' linspace(cfg.col_map{iC}(2),cfg.col_map{iC+1}(2),ceil(1000/(numel(cfg.col_map)-1)))' linspace(cfg.col_map{iC}(3),cfg.col_map{iC+1}(3),ceil(1000/(numel(cfg.col_map)-1)))']; ];
        end
        
        % Scale scores
        neu_dta_use = cell2mat( neu_dev_dta( grp.tle_post_3T_ATLonly_left, neu_dev_col ) );
        xdt_dta_use = cell2mat( cog_dta( grp.tle_post_3T_ATLonly_left, cog_col_one ));
        ydt_dta_use = cell2mat( cog_dta( grp.tle_post_3T_ATLonly_left, cog_col_two ));
        min_val = min(neu_dta_use);
        max_val = max(neu_dta_use);
        min_val = min_val - ((max_val-min_val)*.05);
        
        % Make input data
        clear xdt ydt fce_col
        for iCL = 1:numel(neu_dta_use)
            xdt{iCL} = xdt_dta_use(iCL);
            ydt{iCL} = ydt_dta_use(iCL);
            if ~isnan(neu_dta_use(iCL))
                fce_col{iCL} = col_map(ceil(((neu_dta_use(iCL)-min_val) / (max_val - min_val))*1000)-1,:);
            else
                fce_col{iCL} = rgb('white');
            end
        end
        
        % Plot
        fcfg = [];
        
        fcfg.ydt     = xdt;
        fcfg.xdt     = ydt;
        
        fcfg.fce_col = fce_col;
        fcfg.edg_col = repmat({rgb('black')},1,numel(xdt));
        
        fcfg.ylb = { mmil_spec_char(cog_dta_one{1},{'_'},{' '}) };
        fcfg.xlb = { mmil_spec_char(cog_dta_two{1},{'_'},{' '})  };
        
        fcfg.hln     = -1.28;
        fcfg.hln_col = rgb('red');
        
        fcfg.vln     = -1.28;
        fcfg.vln_col = rgb('red');
        
        fcfg.jtr = 0;
        fcfg.ylm = [ -5 1 ];
        fcfg.xlm = [ -5 1 ];
        
        fcfg.out_dir = [ out_dir '/' 'phenotype_scatters_dev' '/' ];
        fcfg.out_nme = [ roi_sub_mse{iD} '_' neu_238_dta_col{neu_238_col} '_BY_' cog_dta_one{1} ];
        
        ejk_scatter(fcfg)
        
        % Post-operative Phenotype Scatter (LTLE) - 238 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        % Scale scores
        neu_dta_use = cell2mat( neu_238_dta( grp.tle_post_3T_ATLonly_left, neu_238_col ) );
        xdt_dta_use = cell2mat( cog_dta( grp.tle_post_3T_ATLonly_left, cog_col_one ));
        ydt_dta_use = cell2mat( cog_dta( grp.tle_post_3T_ATLonly_left, cog_col_two ));
        min_val = min(neu_dta_use);
        max_val = max(neu_dta_use);
        min_val = min_val - ((max_val-min_val)*.05);
        
        % Make input data
        clear xdt ydt fce_col
        for iCL = 1:numel(neu_dta_use)
            xdt{iCL} = xdt_dta_use(iCL);
            ydt{iCL} = ydt_dta_use(iCL);
            if ~isnan(neu_dta_use(iCL))
                fce_col{iCL} = col_map(ceil(((neu_dta_use(iCL)-min_val) / (max_val - min_val))*1000)-1,:);
            else
                fce_col{iCL} = rgb('white');
            end
        end
        
        % Plot
        fcfg = [];
        
        fcfg.ydt     = xdt;
        fcfg.xdt     = ydt;
        
        fcfg.fce_col = fce_col;
        fcfg.edg_col = repmat({rgb('black')},1,numel(xdt));
        
        fcfg.ylb = { mmil_spec_char(cog_dta_one{1},{'_'},{' '}) };
        fcfg.xlb = { mmil_spec_char(cog_dta_two{1},{'_'},{' '})  };
        
        fcfg.hln     = -1.28;
        fcfg.hln_col = rgb('red');
        
        fcfg.vln     = -1.28;
        fcfg.vln_col = rgb('red');
        
        fcfg.jtr = 0;
        fcfg.ylm = [ -5 1 ];
        fcfg.xlm = [ -5 1 ];
        
        fcfg.out_dir = [ out_dir '/' 'phenotype_scatters_238' '/' ];
        fcfg.out_nme = [ roi_sub_mse{iD} '_' neu_238_dta_col{neu_238_col} '_BY_' cog_dta_one{1} ];
        
        ejk_scatter(fcfg)
        
    end
end

%% Whole Group - DEV focus
ejk_chk_dir([ out_dir '/' 'dev_only' '/' 'post_scatters' '/' ])
ejk_chk_dir([ out_dir '/' 'dev_only' '/' 'pre_scatters' '/' ])
ejk_chk_dir([ out_dir '/' 'dev_only' '/' 'phenotype_scatters' '/' ])
cog_dta_one = { 'log_mem_nor_scr_two_pst' }; % 
cog_dta_two = { 'log_mem_nor_scr_two_zsc' }; % 

roi_sub_fle = { 'subcort_vol_XREP_norm_IntracranialVolume_QC.csv' ...
                'subcort_vol_XREP_LateralityIndex_QC.csv' ...
                'fiber_FA_XREP_QC.csv' ...
                'fiber_FA_XREP_LateralityIndex_QC.csv' ...
                'wmparc_FA_wm_aparc_annot_XREP_QC.csv' ...
                'wmparc_FA_wm_aparc_annot_XREP_LateralityIndex_QC.csv' };
roi_sub_mse = { 'volume'       ...
                'volume_LI'    ...
                'fiber_FA'     ...
                'fiber_FA_LI'  ...
                'wmparc_FA_wm' ...
                'wmparc_FA_wm_LI' };
roi_sub_nme = { { 'Left_Hippocampus' } ...
                { 'Hippocampus' } ...
                { 'L_Unc' 'L_ILF' } ...
                { 'Unc' 'ILF' } ...                
                { 'lh_entorhinal' 'lh_lateralorbitofrontal' 'lh_medialorbitofrontal' 'lh_parsopercularis' } ...
                { 'entorhinal' 'lateralorbitofrontal' 'medialorbitofrontal' 'parsopercularis'  } };
  
fcfg = [];
fcfg.dta_loc = [ dta_fld '/' 'Cognitive_QC.csv' ];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

for iD = 1:numel(roi_sub_fle)

    fcfg = [];
    fcfg.dta_loc = [dta_fld '/' strrep(roi_sub_fle{iD},'XREP','dev') ];
    fcfg.dta_col = 5;
    [ neu_dev_dta, neu_dev_dta_sbj, neu_dev_dta_col] = ejk_dta_frm( fcfg );    
    
    for iR = 1:numel(roi_sub_nme{iD})
        
        neu_dev_col = find(strcmpi(neu_dev_dta_col,roi_sub_nme{iD}{iR}));
               
        % Pre-operative (all)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cog_col_two = find(strcmpi(cog_dta_col,cog_dta_two{1}));
        
        fcfg = [];
        
        fcfg.xdt     = { cell2mat(neu_dev_dta( grp.tle_controls_pre_3T_allSurg_all , neu_dev_col )) };
        fcfg.ydt     = { cell2mat(cog_dta( grp.tle_controls_pre_3T_allSurg_all, cog_col_two )) };
        
        fcfg.fce_col = { rgb('burnt orange') };
        fcfg.edg_col = { rgb('black') };
        
        fcfg.trd_lne = [ 1 1 ];
        
        fcfg.ylb = { mmil_spec_char(cog_dta_two{1},{'_'},{' '}) };
        fcfg.xlb = { mmil_spec_char(neu_dev_dta_col{neu_dev_col},{'_'},{' '})  };
        
        fcfg.hln     = -1.28;
        fcfg.hln_col = rgb('red');
        
        fcfg.vln     = -1.28;
        fcfg.vln_col = rgb('red');
        
        fcfg.ylm = [ -5 1 ];
        
        fcfg.out_dir = [ out_dir '/' 'dev_only' '/' 'pre_scatters' '/' ];
        fcfg.out_nme = [ roi_sub_mse{iD} '_' neu_dev_dta_col{neu_dev_col} '_BY_' cog_dta_two{1} ];
        
        ejk_scatter(fcfg)        
        
        % Post-operative scatter (LTLE)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cog_col_one = find(strcmpi(cog_dta_col,cog_dta_one{1}));
        
        fcfg = [];
        
        fcfg.xdt     = { cell2mat(neu_dev_dta( grp.tle_post_3T_ATLonly_left , neu_dev_col )) };
        fcfg.ydt     = { cell2mat(cog_dta( grp.tle_post_3T_ATLonly_left, cog_col_one )) };
        
        fcfg.fce_col = { rgb('burnt orange') };
        fcfg.edg_col = { rgb('black') };
        
        fcfg.trd_lne = [ 1 1 ];
        
        fcfg.ylb = { mmil_spec_char(cog_dta_one{1},{'_'},{' '}) };
        fcfg.xlb = { mmil_spec_char(neu_dev_dta_col{neu_dev_col},{'_'},{' '})  };
        
        fcfg.hln     = -1.28;
        fcfg.hln_col = rgb('red');
        
        fcfg.vln     = -1.28;
        fcfg.vln_col = rgb('red');
        
        fcfg.ylm = [ -5 1 ];
        
        fcfg.out_dir = [ out_dir '/' 'dev_only' '/' 'post_scatters' '/' ];
        fcfg.out_nme = [ roi_sub_mse{iD} '_' neu_dev_dta_col{neu_dev_col} '_BY_' cog_dta_one{1} ];
        
        ejk_scatter(fcfg)        
        
        % Post-operative Phenotype Scatter (LTLE) - DEV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make color map
        cfg.col_map = { rgb('neon red') rgb('red') rgb('grey')-0.2 rgb('blue') rgb('neon blue') };
        col_map = [];
        for iC = 1:numel(cfg.col_map)-1
            col_map = [col_map ; [linspace(cfg.col_map{iC}(1),cfg.col_map{iC+1}(1),ceil(1000/(numel(cfg.col_map)-1)))' linspace(cfg.col_map{iC}(2),cfg.col_map{iC+1}(2),ceil(1000/(numel(cfg.col_map)-1)))' linspace(cfg.col_map{iC}(3),cfg.col_map{iC+1}(3),ceil(1000/(numel(cfg.col_map)-1)))']; ];
        end
        
        % Scale scores
        neu_dta_use = cell2mat( neu_dev_dta( grp.tle_post_3T_ATLonly_left, neu_dev_col ) );
        xdt_dta_use = cell2mat( cog_dta( grp.tle_post_3T_ATLonly_left, cog_col_one ));
        ydt_dta_use = cell2mat( cog_dta( grp.tle_post_3T_ATLonly_left, cog_col_two ));
        min_val = min(neu_dta_use);
        max_val = max(neu_dta_use);
        min_val = min_val - ((max_val-min_val)*.05);
        
        % Make input data
        clear xdt ydt fce_col
        [ ~, srt_ind ] = sort(neu_dta_use);
        for iCL = 1:numel(neu_dta_use)
            xdt{iCL} = xdt_dta_use(iCL);
            ydt{iCL} = ydt_dta_use(iCL);
            if ~isnan(neu_dta_use(iCL))
                fce_col{iCL} = col_map(floor((find(srt_ind==iCL) / numel(neu_dta_use)*1000)),:);% col_map(ceil(((neu_dta_use(iCL)-min_val) / (max_val - min_val))*1000)-1,:);
            else
                fce_col{iCL} = rgb('white');
            end
        end
        
        % Plot
        fcfg = [];
        
        fcfg.ydt     = xdt;
        fcfg.xdt     = ydt;
        
        fcfg.fce_col = fce_col;
        fcfg.edg_col = repmat({rgb('black')},1,numel(xdt));
        
        fcfg.ylb = { mmil_spec_char(cog_dta_one{1},{'_'},{' '}) };
        fcfg.xlb = { mmil_spec_char(cog_dta_two{1},{'_'},{' '})  };
        
        fcfg.hln     = -1.28;
        fcfg.hln_col = rgb('red');
        
        fcfg.vln     = -1.28;
        fcfg.vln_col = rgb('red');
        
        fcfg.jtr = 0;
        fcfg.ylm = [ -5 1 ];
        fcfg.xlm = [ -5 1 ];
        
        fcfg.out_dir = [ out_dir '/' 'dev_only' '/' 'phenotype_scatters' '/' ];
        fcfg.out_nme = [ roi_sub_mse{iD} '_' neu_dev_dta_col{neu_dev_col} '_BY_' cog_dta_one{1} ];
        
        ejk_scatter(fcfg)
                
    end
end

















            


            


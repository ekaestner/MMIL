idv_out_dir = [ out_dir '/' 'individual_revision' '/']; ejk_chk_dir(idv_out_dir);

cor_nme = { 'Age' 'Age at Onset' 'Duration' 'Education Years' };
grp_nme = { 'Gender' 'Group' 'Handedness' 'Side of Epilepsy' 'Site' 'MTS Status' };

%%
load([ dta_dir '/' 'demographics.mat' ]);
load([ dta_dir '/' 'prediction.mat' ]);

%% Total performance


%% Create grp variable
grp.main.HC    = find(strcmpi(cmb_dta(:,strcmpi(cmb_col,'Group')),'HC'));
grp.main.TLE   = find(strcmpi(cmb_dta(:,strcmpi(cmb_col,'Group')),'TLE'));
grp.main.MRI_NEG    = find(strcmpi(cmb_dta(:,strcmpi(cmb_col,'MTS Status')),'MRI-'));
grp.main.MRI_POS   = find(strcmpi(cmb_dta(:,strcmpi(cmb_col,'MTS Status')),'MRI+'));
grp.main.TOTAL = [ grp.main.HC ; grp.main.TLE ];

grp_man_nme = fieldnames(grp.main);

grp_row_nme     = { 'r-value' 'p-value' 'test' };
grp_row_nme_hld = [];
for iG = 1:numel(grp_man_nme)
    grp_row_nme_hld = [ grp_row_nme_hld strcat(grp_man_nme{iG},'_',grp_row_nme) ];
end

%%
rho_tst_tbl = cell(numel([cor_nme grp_nme ]),3*numel(grp_man_nme));
rvl_tst_tbl = cell(numel([cor_nme grp_nme ]),3*numel(grp_man_nme));

%% Correlations: Age, Education, Age of Onset, Duration
col_ind = 1;
for iG = 1:numel(grp_man_nme)

    tst_ind = 1;
    
    % Collate
    sbj_use_ind = zeros(size(cmb_dta,1),1);
    sbj_use_ind(grp.main.(grp_man_nme{iG})) = 1;
    
    % Split by patient population
    for iCR = 1:numel(cor_nme)
        
        % Collate        
        cor_col = strcmpi(cmb_col,cor_nme{iCR});
        cor_row = ~cellfun( @isempty, cmb_dta(:,cor_col)) & sbj_use_ind;
        
        cor_dta     = cell2mat(cmb_dta(cor_row,cor_col));
        cor_sbj     = cmb_sbj(cor_row,1);
        prd_scr_use = prd_scr(cor_row,1);
        
        [ rho, poh ] = corr(cor_dta,prd_scr_use,'Type','Spearman','rows','complete');
        [ rvl, pvl ] = corr(cor_dta,prd_scr_use,'Type','Pearson','rows','complete');
        
        % Plot
        fcfg = [];
        
        fcfg.xdt = { cor_dta };
        fcfg.ydt = { prd_scr_use };
        
        fcfg.fce_col     = { rgb('light orange') };
        fcfg.edg_col     = { [0 0 0]              };
        fcfg.box_plt_col = { rgb('dark orange')  };
        
        fcfg.box_plt = ones(1,numel(fcfg.xdt));
        fcfg.xlb = { cor_nme{iCR} };
        fcfg.ylb = {'Accuracy'};
        fcfg.ylm = [ 0 100 ];
        
        fcfg.ttl = [ cor_nme{iCR} ' ' 'pvl:' ' ' num2str(roundsd(pvl,2)) ' ' num2str(roundsd(poh,2)) ];
        
        fcfg.trd_lne = [1 1];
        
        fcfg.out_dir = idv_out_dir;
        fcfg.out_nme = [ 'corr' '_' 'p' num2str(iCR) '_' mmil_spec_char(cor_nme{iCR},{' '},{'_'}) '_' grp_man_nme{iG} ];
        
        ejk_scatter(fcfg)
        
        %
        rho_str = num2str(roundsd(rho,2));
        poh_str = num2str(roundsd(poh,2));
        if strcmpi(rho_str(1),'-')
            rho_tst_tbl(tst_ind,col_ind:col_ind+2) = [ {rho} {poh} [ 'r=' rho_str([1 3:end]) '; p=' poh_str(2:end)] ];
        else
            rho_tst_tbl(tst_ind,col_ind:col_ind+2) = [ {rho} {poh} [ 'r=' rho_str(2:end) '; p=' poh_str(2:end)] ];
        end
        
        rvl_str = num2str(roundsd(rvl,2));
        pvl_str = num2str(roundsd(pvl,2));
        if strcmpi(rvl_str(1),'-')
            rvl_tst_tbl(tst_ind,col_ind:col_ind+2) = [ {rvl} {pvl} [ 'r=' rvl_str([1 3:end]) '; p=' pvl_str(2:end)] ];
        else
            rvl_tst_tbl(tst_ind,col_ind:col_ind+2) = [ {rvl} {pvl} [ 'r=' rvl_str(2:end) '; p=' pvl_str(2:end)] ];
        end
                
        tst_ind = tst_ind + 1;
        
    end
    col_ind = col_ind + 3;
end

%% Group comparisons: Gender, Side of Epilepsy, Handedness,
hld_tst_ind = tst_ind;
col_ind = 1;
for iG = 1:numel(grp_man_nme)
    
    tst_ind = hld_tst_ind;
    
    % Collate
    sbj_use_ind = zeros(size(cmb_dta,1),1);
    sbj_use_ind(grp.main.(grp_man_nme{iG})) = 1;
    
    for iAN = 1:numel(grp_nme)
        
        % Collate
        anv_col = strcmpi(cmb_col,grp_nme{iAN});
        anv_row = ~strcmpi( cmb_dta(:,anv_col),'')  & sbj_use_ind;
        
        anv_dta     = cmb_dta(anv_row,anv_col);
        cor_sbj     = cmb_sbj(anv_row,1);
        prd_scr_use = prd_scr(anv_row,1);
        
        [ pvl_hld, tst_hld, stt_hld ] = anova1(prd_scr_use, anv_dta,'off');
        
        % Plot
        lvl_nme = unique(anv_dta);
        
        if ~isempty(lvl_nme)
            fcfg = [];
            
            for iL = 1:numel(lvl_nme)
                fcfg.xdt{iL} = iL;
                fcfg.ydt{iL} = prd_scr_use(strcmpi(anv_dta,lvl_nme{iL}));
                fcfg.xlb{iL} = lvl_nme{iL};
            end
            
            fcfg.fce_col     = repmat( {[0 0 0]}, 1, numel(fcfg.xdt) );
            fcfg.edg_col     = repmat( {[0 0 0]}, 1, numel(fcfg.xdt) );
            fcfg.box_plt_col = repmat( {[0.6 0.6 0.6]}, 1, numel(fcfg.xdt) );
            
            fcfg.box_plt = ones(1,numel(fcfg.xdt));
            
            fcfg.ylb = {'Accuracy'};
            fcfg.ylm = [ 0 100 ];
            
            fcfg.mkr_sze = repmat(15,1,numel(fcfg.xdt));
            fcfg.aph_val = 0.45;
            
            fcfg.ttl = [ grp_nme{iAN} ' ' 'pvl:' ' ' num2str(roundsd(pvl_hld,2)) ];
            
            fcfg.out_dir = idv_out_dir;
            fcfg.out_nme = [ 'anova' '_' 'p' num2str(iAN) '_' mmil_spec_char(grp_nme{iAN},{' '},{'_'}) '_' grp_man_nme{iG} ];
            
            ejk_scatter(fcfg)
            
            %
            fvl_str = num2str(roundsd(tst_hld{2,5},3));
            pvl_str = num2str(roundsd(pvl_hld,2));
            rho_tst_tbl(tst_ind,col_ind:col_ind+2) = [ tst_hld(2,5) {pvl_hld} [ 'F(' num2str(stt_hld.df) ')=' fvl_str '; p=' pvl_str(2:end)]];
            rvl_tst_tbl(tst_ind,col_ind:col_ind+2) = [ tst_hld(2,5) {pvl_hld} [ 'F(' num2str(stt_hld.df) ')=' fvl_str '; p=' pvl_str(2:end)]];
        end
        
        tst_ind = tst_ind + 1;
    end
    col_ind = col_ind + 3;
end

%%
cell2csv( [ idv_out_dir '/' 'rho_tbl.csv' ], [ {''} grp_row_nme_hld ; [cor_nme' ; grp_nme' ] rho_tst_tbl ] )
cell2csv( [ idv_out_dir '/' 'rvl_tbl.csv' ], [ {''} grp_row_nme_hld ; [cor_nme' ; grp_nme' ] rvl_tst_tbl ] )
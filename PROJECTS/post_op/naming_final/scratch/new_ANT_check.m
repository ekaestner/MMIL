clear; clc;

prj_dir = '/home/ekaestne/PROJECTS/OUTPUT';

new_prj_nme = 'PostOperative/Naming_final_sample';
org_prj_nme = 'PostOperative/Naming_final';

new_grp = load( [ prj_dir '/' new_prj_nme '/' 'groups.mat' ] );
old_grp = load( [ prj_dir '/' org_prj_nme '/' 'groups.mat' ] );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' new_prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
[ new_cog_dta, new_cog_dta_sbj, new_cog_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' org_prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
[ old_cog_dta, old_cog_dta_sbj, old_cog_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' org_prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' org_prj_nme '/' 'Data' '/' 'wmparc_FA_wm_aparc_annot_QC.csv'];
fcfg.dta_col = 2;
[ wmp_dta, wmp_sbj, wmp_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' org_prj_nme '/' 'Data' '/' 'fiber_FA_QC.csv'];
fcfg.dta_col = 2;
[ trc_dta, trc_sbj, trc_col] = ejk_dta_frm( fcfg );

neu_dta = [ wmp_dta trc_dta ];
neu_col = [ wmp_col trc_col ];
neu_col_ind(1) = find(strcmpi(neu_col,'rh_fusiform'));
neu_col_ind(2) = find(strcmpi(neu_col,'L_ILF'));
neu_col_ind(3) = find(strcmpi(neu_col,'L_IFO'));
neu_ylm = {[.25 .4] [.35 .6] [.35 .6]};

cln_dta = [ cln_dta new_cog_dta ];
cln_col = [ cln_dta_col new_cog_dta_col ];
cln_col_ind = [ 16 4 5 ];
cln_ylm = {[35 51] [15 65] [8 20]};

%% Correlation of ANT & BNT
grp_nme = { 'tle_post_3T_ATLonly_all' 'tle_post_3T_ATLonly_left' 'tle_post_3T_ATLonly_right' };

figure()

for iG = 1:numel(grp_nme)

    subplot(3,3,iG)
    
    xdt = cell2mat(new_cog_dta(new_grp.grp.(grp_nme{iG}),3));
    ydt = cell2mat(new_cog_dta(new_grp.grp.(grp_nme{iG}),4));
    scatter( xdt, ydt, 'r', 'filled'); hold on;
    use_idx = ~isnan(xdt) & ~isnan(ydt);
    trd_lne_fit = polyfit(xdt(use_idx),ydt(use_idx),1);
    trd_lne_fit = polyval(trd_lne_fit,xdt(use_idx));
    plot( xdt(use_idx), trd_lne_fit, 'Color', rgb('red'),'LineWidth',2)  
    
    xdt = cell2mat(old_cog_dta(old_grp.grp.(grp_nme{iG}),3));
    ydt = cell2mat(old_cog_dta(old_grp.grp.(grp_nme{iG}),4));
    scatter( xdt, ydt, 'b', 'filled'); hold on;
    use_idx = ~isnan(xdt) & ~isnan(ydt);
    trd_lne_fit = polyfit(xdt(use_idx),ydt(use_idx),1);
    trd_lne_fit = polyval(trd_lne_fit,xdt(use_idx));
    plot( xdt(use_idx), trd_lne_fit, 'Color', rgb('blue'),'LineWidth',2)
    
    title( mmil_spec_char(grp_nme{iG},{'_'},{' '}) ); xlabel('BNT'); ylabel('ANT'); ylim([-5 5]); xlim([-5 5]);
    
end

tightfig();

%% Neurobio
grp_nme = { 'tle_post_3T_ATLonly_all' 'tle_post_3T_ATLonly_left' 'tle_post_3T_ATLonly_right' };

% R-Fusiform, L-ILF, L-IFOF x All/Left/Right x ANT
fig_ind = 1;
figure()
for iG = 1:numel(grp_nme)
    for iC = 1:numel(neu_col_ind)
        
        subplot(3,3,fig_ind)
        
        xdt = cell2mat(new_cog_dta(new_grp.grp.(grp_nme{iG}),4));
        ydt = cell2mat(neu_dta(new_grp.grp.(grp_nme{iG}),neu_col_ind(iC)));
        scatter( xdt, ydt, 'r', 'filled'); hold on;
        use_idx = ~isnan(xdt) & ~isnan(ydt);
        trd_lne_fit = polyfit(xdt(use_idx),ydt(use_idx),1);
        trd_lne_fit = polyval(trd_lne_fit,xdt(use_idx));
        plot( xdt(use_idx), trd_lne_fit, 'Color', rgb('red'),'LineWidth',2)
        
        xdt = cell2mat(old_cog_dta(old_grp.grp.(grp_nme{iG}),4));
        ydt = cell2mat(neu_dta(old_grp.grp.(grp_nme{iG}),neu_col_ind(iC)));
        scatter( xdt, ydt, 'b', 'filled'); hold on;
        use_idx = ~isnan(xdt) & ~isnan(ydt);
        trd_lne_fit = polyfit(xdt(use_idx),ydt(use_idx),1);
        trd_lne_fit = polyval(trd_lne_fit,xdt(use_idx));
        plot( xdt(use_idx), trd_lne_fit, 'Color', rgb('blue'),'LineWidth',2)
        
        title( mmil_spec_char(grp_nme{iG},{'_'},{' '}) ); xlabel('BNT'); ylabel( mmil_spec_char(neu_col{neu_col_ind(iC)},{'_'},{' '})); ylim(neu_ylm{iC}); xlim([-5 5]);
        
        fig_ind = fig_ind + 1;
        
    end
end

tightfig();

%% Clinical
grp_nme = { 'tle_post_3T_ATLonly_all' 'tle_post_3T_ATLonly_left' 'tle_post_3T_ATLonly_right' };

% R-Fusiform, L-ILF, L-IFOF x All/Left/Right x ANT
fig_ind = 1;
figure()
for iG = 1:numel(grp_nme)
    for iC = 1:numel(cln_col_ind)
        
        subplot(3,3,fig_ind)
        
        xdt = cell2mat(new_cog_dta(new_grp.grp.(grp_nme{iG}),4));
        ydt = cell2mat(cln_dta(new_grp.grp.(grp_nme{iG}),cln_col_ind(iC)));
        scatter( xdt, ydt, 'r', 'filled'); hold on;
        use_idx = ~isnan(xdt) & ~isnan(ydt);
        trd_lne_fit = polyfit(xdt(use_idx),ydt(use_idx),1);
        trd_lne_fit = polyval(trd_lne_fit,xdt(use_idx));
        plot( xdt(use_idx), trd_lne_fit, 'Color', rgb('red'),'LineWidth',2)
        
        xdt = cell2mat(old_cog_dta(old_grp.grp.(grp_nme{iG}),4));
        ydt = cell2mat(cln_dta(old_grp.grp.(grp_nme{iG}),cln_col_ind(iC)));
        scatter( xdt, ydt, 'b', 'filled'); hold on;
        use_idx = ~isnan(xdt) & ~isnan(ydt);
        trd_lne_fit = polyfit(xdt(use_idx),ydt(use_idx),1);
        trd_lne_fit = polyval(trd_lne_fit,xdt(use_idx));
        plot( xdt(use_idx), trd_lne_fit, 'Color', rgb('blue'),'LineWidth',2)
        
        title( mmil_spec_char(grp_nme{iG},{'_'},{' '}) ); xlabel('BNT'); ylabel( mmil_spec_char(cln_col{cln_col_ind(iC)},{'_'},{' '})); ylim(cln_ylm{iC}); xlim([-5 5]);
        
        fig_ind = fig_ind + 1;
        
    end
end

tightfig();

%%

[ trc_sbj(new_grp.grp.(grp_nme{2})) new_cog_dta(new_grp.grp.(grp_nme{2}),3) new_cog_dta(new_grp.grp.(grp_nme{2}),4) neu_dta(new_grp.grp.(grp_nme{2}),neu_col_ind(1)) neu_dta(new_grp.grp.(grp_nme{2}),neu_col_ind(2)) ]


[ trc_sbj(new_grp.grp.(grp_nme{3})) new_cog_dta(new_grp.grp.(grp_nme{3}),3) new_cog_dta(new_grp.grp.(grp_nme{3}),4) neu_dta(new_grp.grp.(grp_nme{3}),neu_col_ind(1)) neu_dta(new_grp.grp.(grp_nme{3}),neu_col_ind(2)) ]






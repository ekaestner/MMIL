%% Performance on 2D ComBat
tbl_out_dir = [ out_dir '/' 'Tables' '_' 'Revision' '/' ]; ejk_chk_dir(tbl_out_dir);

load([ dta_dir '/' 'performance_flip_reihaneh_revision.mat' ])

mdl_int_nme = { 'fcn_3dm_org'                    'fcn_3dm_cbt'                'fcn_2dm_org'                    'fcn_2dm_cbt'  };
mes_int_nme = { 'Accuracy' 'Specificity' 'Sensitivity' 'PPV' 'NPV' 'F1'};

avg_out_tbl = cell(numel(mdl_int_nme),numel(mes_int_nme));
med_out_tbl = cell(numel(mdl_int_nme),numel(mes_int_nme));
for iC = 1:numel(mes_int_nme)
    for iR = 1:numel(mdl_int_nme)
        med_val = nanmedian(cell2mat(dta_hld.(mdl_int_nme{iR}).(mes_int_nme{iC})(:,3)));
        avg_val = nanmean(cell2mat(dta_hld.(mdl_int_nme{iR}).(mes_int_nme{iC})(:,3)));
        std_val = nanstd(cell2mat(dta_hld.(mdl_int_nme{iR}).(mes_int_nme{iC})(:,3)));
        iqr_val = iqr(cell2mat(dta_hld.(mdl_int_nme{iR}).(mes_int_nme{iC})(:,3)));
        max_val = nanmax(cell2mat(dta_hld.(mdl_int_nme{iR}).(mes_int_nme{iC})(:,3)));
        min_val = nanmin(cell2mat(dta_hld.(mdl_int_nme{iR}).(mes_int_nme{iC})(:,3)));
        avg_out_tbl{iR,iC} = sprintf('%.1f (%.1f); %.1f',avg_val,std_val, max_val);
        med_out_tbl{iR,iC} = sprintf('%.1f (%.1f); %.1f',med_val,iqr_val, max_val);
    end
end

cell2csv([ tbl_out_dir '/' 'PerformanceTable_mean.csv'],   [ {''}  mes_int_nme; mdl_int_nme' avg_out_tbl ])
cell2csv([ tbl_out_dir '/' 'PerformanceTable_median.csv'], [ {''}  mes_int_nme; mdl_int_nme' med_out_tbl ])

%% Calculate performance on MRI-negative patients
load([ dta_dir '/' 'prediction.mat' ])
load([ dta_dir '/' 'demographics.mat' ])

prd_dir = [ dta_dir '/' 'predictions' '/'];
mdl_nme = { 'predictions_FCNet2D'  }; %'predictions_har_FCNet2D' 'predictions_har_FCNet3D'
run_num = { 'r0' 'r1' 'r2' 'r3' 'r4' 'r5' 'r6' 'r7' 'r8' 'r9' }; % 
fld_num = { 's0' 's1' 's2' 's3' 's4' };

clear prd_dta
prd_dta.sbj_nme = prd_dta_sbj;
prd_dta.predictions_FCNet3D = prd_scr;

grp.main.MRI_NEG    = find(strcmpi(cmb_dta(:,strcmpi(cmb_col,'MTS Status')),'MRI-'));
grp.main.MRI_POS   = find(strcmpi(cmb_dta(:,strcmpi(cmb_col,'MTS Status')),'MRI+'));

for iM = 1:numel(mdl_nme)
    prd_dta.(mdl_nme{iM}) = zeros(numel(prd_dta.sbj_nme),1);
    
    for iR = 1:numel(run_num)
        for iF = 1:numel(fld_num)
            fcfg = [];
            fcfg.dta_loc = [ prd_dir '/' mdl_nme{iM} '/' 'pred' '_' run_num{iR} '_' fld_num{iF} '.csv'];
            [fld_dta, fld_sbj, fld_col] = ejk_dta_frm( fcfg );
             
            for iS = 1:numel(fld_sbj)
                prd_dta.(mdl_nme{iM})( strcmpi(prd_dta.sbj_nme, fld_sbj{iS}) ) = prd_dta.(mdl_nme{iM})( strcmpi(prd_dta.sbj_nme, fld_sbj{iS}) ) + fld_dta{iS,1};
            end            
        end
    end
    
    prd_dta.(mdl_nme{iM}) = (prd_dta.(mdl_nme{iM}) / 10) * 100;    
end

ttt = [ mean(prd_dta.predictions_FCNet3D(grp.main.MRI_POS)) mean(prd_dta.predictions_FCNet3D(grp.main.MRI_NEG)) ; ...
  mean(prd_dta.predictions_FCNet2D(grp.main.MRI_POS)) mean(prd_dta.predictions_FCNet2D(grp.main.MRI_NEG)) ]; %; ... 
  %mean(prd_dta.predictions_har_FCNet3D(grp.main.MRI_POS)) mean(prd_dta.predictions_har_FCNet3D(grp.main.MRI_NEG)) ; ...
  %mean(prd_dta.predictions_har_FCNet2D(grp.main.MRI_POS)) mean(prd_dta.predictions_har_FCNet2D(grp.main.MRI_NEG)) ]
ttt = [ ttt ; [ttt(1,1) - ttt(2,1)] [ttt(1,2) - ttt(2,2)] ]; 
ttt

%% Correlation on MRI-negative patients
cnn_3dm_mri_neg_revision

%% Alternate model comparison test
load([ dta_dir '/' 'performance_flip_reihaneh_revision.mat' ])

% 3D_org vs 2D_org Wilcoxon
signrank( cell2mat(dta_hld.fcn_3dm_org.Accuracy(:,3)), cell2mat(dta_hld.fcn_2dm_org.Accuracy(:,3)) )

% 3D_org vs 3D_cbt Wilcoxon
signrank( cell2mat(dta_hld.fcn_3dm_org.Accuracy(:,3)), cell2mat(dta_hld.fcn_3dm_cbt.Accuracy(:,3)) )

% 3D_cbt vs 2D_cbt Wilcoxon
signrank( cell2mat(dta_hld.fcn_3dm_cbt.Accuracy(:,3)), cell2mat(dta_hld.fcn_2dm_cbt.Accuracy(:,3)) )

sum( cell2mat(dta_hld.fcn_3dm_cbt.Accuracy(:,3)) > cell2mat(dta_hld.fcn_2dm_cbt.Accuracy(:,3)) )
(sum( cell2mat(dta_hld.fcn_3dm_cbt.Accuracy(:,3)) > cell2mat(dta_hld.fcn_2dm_cbt.Accuracy(:,3)) ) / 50) * 100

% 2D_org vs 2D_cbt Wilcoxon
signrank( cell2mat(dta_hld.fcn_2dm_org.Accuracy(:,3)), cell2mat(dta_hld.fcn_2dm_cbt.Accuracy(:,3)) )

sum( cell2mat(dta_hld.fcn_2dm_org.Accuracy(:,3)) > cell2mat(dta_hld.fcn_2dm_cbt.Accuracy(:,3)) )
(sum( cell2mat(dta_hld.fcn_2dm_org.Accuracy(:,3)) > cell2mat(dta_hld.fcn_2dm_cbt.Accuracy(:,3)) ) / 50) * 100


% signrank( cell2mat(dta_hld.fcn_3dm_org.Accuracy(:,3)), cell2mat(dta_hld.fcn_2dm_cbt.Accuracy(:,3)) )
% signrank( cell2mat(dta_hld.fcn_2dm_cbt.Accuracy(:,3)), cell2mat(dta_hld.fcn_3dm_org.Accuracy(:,3)) )
% signrank( cell2mat(dta_hld.fcn_2dm_cbt.Accuracy(:,3)), cell2mat(dta_hld.fcn_2dm_org.Accuracy(:,3)) )



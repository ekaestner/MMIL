clear; clc;

%% DTI
mps_238 = mmil_readtext('/home/mmilmcdRSI/MetaData/MCD_RSI/ROI_Summaries/DTI_all.csv');
mps_dev = mmil_readtext('/space/mcdonald-syn01/1/data/MCD_MRI/MetaData/MCD_MRI/ROI_Summaries/DTI_all.csv');

mps_238_mri = mmil_readtext('/home/mmilmcdRSI/MetaData/MCD_RSI/ROI_Summaries/MRI_all.csv');
mps_dev_mri = mmil_readtext('/space/mcdonald-syn01/1/data/MCD_MRI/MetaData/MCD_MRI/ROI_Summaries/MRI_all.csv');

%%
[ ovr_sbj, ovr_sbj_ind_238, ovr_sbj_ind_dev ] = intersect( mps_238(:,2), mps_dev(:,2) );
ovr_sbj(1:3,:)=[];ovr_sbj_ind_238(1:3)=[];ovr_sbj_ind_dev(1:3)=[];

[ ovr_sbj mps_238(ovr_sbj_ind_238,2) mps_dev(ovr_sbj_ind_dev,2)];

wmp_wfa_col_238 = string_find(mps_238(1,:),{'wmparc_FA-wm'});
wmp_wfa_col_dev = string_find(mps_dev(1,:),{'wmparc_FA-wm'});

wmp_wmd_col_238 = string_find(mps_238(1,:),{'wmparc_MD-wm'});
wmp_wmd_col_dev = string_find(mps_dev(1,:),{'wmparc_MD-wm'});

trc_tfa_col_238 = string_find(mps_238(1,:),{'fiber_FA'});
trc_tfa_col_dev = string_find(mps_dev(1,:),{'fiber_FA'});

trc_tmd_col_238 = string_find(mps_238(1,:),{'fiber_MD'});
trc_tmd_col_dev = string_find(mps_dev(1,:),{'fiber_MD'});

mri_thk_col_238 = string_find(mps_238_mri(1,:),{'thick-ctx'});
    mri_thk_col_238(string_find(mps_238_mri(1,mri_thk_col_238),{'fuzzy'})) = [];
mri_thk_col_dev = string_find(mps_dev_mri(1,:),{'thick-ctx'});
    mri_thk_col_dev(string_find(mps_dev_mri(1,mri_thk_col_dev),{'fuzzy'})) = [];

mri_vol_col_238 = string_find(mps_238(1,:),{''});
mri_vol_col_dev = string_find(mps_dev(1,:),{''});

%% WMPARC FA
ejk_chk_dir('/home/ekaestne/PROJECTS/OUTPUT/mmps_compare/wmparc_FA/figs')
ejk_chk_dir('/home/ekaestne/PROJECTS/OUTPUT/mmps_compare/wmparc_FA/diffs')

for iTF = 1:numel(wmp_wfa_col_238)
    
    % Figure
    figure('Visible','off')
    scatter( cell2mat(mps_238(ovr_sbj_ind_238,wmp_wfa_col_238(iTF))), cell2mat(mps_dev(ovr_sbj_ind_dev,wmp_wfa_col_dev(iTF))), 'filled' )
    title( mmil_spec_char(mps_238{1,wmp_wfa_col_238(iTF)},{'_' '-'},{' ' ' '}) )
    xlabel('MMPS 238')
    ylabel('MMPS DEV')
    print(gcf,[ '/home/ekaestne/PROJECTS/OUTPUT/mmps_compare/wmparc_FA/figs' '/' mps_dev{1,wmp_wfa_col_dev(iTF)} '_' mps_238{1,wmp_wfa_col_238(iTF) } '.png' ],'-dpng')
    close all
    
    % Diffs
    clear avg_dff
    for iS = 1:numel(ovr_sbj_ind_238)
        avg_dff{iS,1} = mps_238{ovr_sbj_ind_238(iS),2};
        avg_dff{iS,2} = mps_dev{ovr_sbj_ind_dev(iS),2};
        avg_dff{iS,3} = nanmean( cell2mat(mps_238(ovr_sbj_ind_238(iS),wmp_wfa_col_238(iTF))) - ...
                                 cell2mat(mps_dev(ovr_sbj_ind_dev(iS),wmp_wfa_col_dev(iTF))) );
    end
    cell2csv([ '/home/ekaestne/PROJECTS/OUTPUT/mmps_compare/wmparc_FA/diffs' '/' mps_dev{1,wmp_wfa_col_dev(iTF)} '_' mps_238{1,wmp_wfa_col_238(iTF) } '.csv' ],avg_dff);
    
end

%% WMPARC MD
ejk_chk_dir('/home/ekaestne/PROJECTS/OUTPUT/mmps_compare/wmparc_MD/figs')
ejk_chk_dir('/home/ekaestne/PROJECTS/OUTPUT/mmps_compare/wmparc_MD/diffs')

for iTF = 1:numel(wmp_wmd_col_238)
    
    % Figure
    figure('Visible','off')
    scatter( cell2mat(mps_238(ovr_sbj_ind_238,wmp_wmd_col_238(iTF))), cell2mat(mps_dev(ovr_sbj_ind_dev,wmp_wmd_col_dev(iTF))), 'filled' )
    title( mmil_spec_char(mps_238{1,wmp_wmd_col_238(iTF)},{'_' '-'},{' ' ' '}) )
    xlabel('MMPS 238')
    ylabel('MMPS DEV')
    print(gcf,[ '/home/ekaestne/PROJECTS/OUTPUT/mmps_compare/wmparc_MD/figs' '/' mps_dev{1,wmp_wmd_col_dev(iTF)} '_' mps_238{1,wmp_wmd_col_238(iTF) } '.png' ],'-dpng')
    close all
    
    % Diffs
    clear avg_dff
    for iS = 1:numel(ovr_sbj_ind_238)
        avg_dff{iS,1} = mps_238{ovr_sbj_ind_238(iS),2};
        avg_dff{iS,2} = mps_dev{ovr_sbj_ind_dev(iS),2};
        avg_dff{iS,3} = nanmean( cell2mat(mps_238(ovr_sbj_ind_238(iS),wmp_wmd_col_238(iTF))) - ...
                                 cell2mat(mps_dev(ovr_sbj_ind_dev(iS),wmp_wmd_col_dev(iTF))) );
    end
    cell2csv([ '/home/ekaestne/PROJECTS/OUTPUT/mmps_compare/wmparc_MD/diffs' '/' mps_dev{1,wmp_wmd_col_dev(iTF)} '_' mps_238{1,wmp_wmd_col_238(iTF) } '.csv' ],avg_dff);
    
end

%% TRACT FA
ejk_chk_dir('/home/ekaestne/PROJECTS/OUTPUT/mmps_compare/tract_FA/figs')
ejk_chk_dir('/home/ekaestne/PROJECTS/OUTPUT/mmps_compare/tract_FA/diffs')

for iTF = 1:numel(trc_tfa_col_238)
    
    % Figure
    figure('Visible','off')
    scatter( cell2mat(mps_238(ovr_sbj_ind_238,trc_tfa_col_238(iTF))), cell2mat(mps_dev(ovr_sbj_ind_dev,trc_tfa_col_dev(iTF))), 'filled' )
    title( mmil_spec_char(mps_238{1,trc_tfa_col_238(iTF)},{'_' '-'},{' ' ' '}) )
    xlabel('MMPS 238')
    ylabel('MMPS DEV')
    print(gcf,[ '/home/ekaestne/PROJECTS/OUTPUT/mmps_compare/tract_FA/figs' '/' mps_dev{1,trc_tfa_col_dev(iTF)} '_' mps_238{1,trc_tfa_col_238(iTF) } '.png' ],'-dpng')
    close all
    
    % Diffs
    clear avg_dff
    for iS = 1:numel(ovr_sbj_ind_238)
        avg_dff{iS,1} = mps_238{ovr_sbj_ind_238(iS),2};
        avg_dff{iS,2} = mps_dev{ovr_sbj_ind_dev(iS),2};
        avg_dff{iS,3} = nanmean( cell2mat(mps_238(ovr_sbj_ind_238(iS),trc_tfa_col_238(iTF))) - ...
                                 cell2mat(mps_dev(ovr_sbj_ind_dev(iS),trc_tfa_col_dev(iTF))) );
    end
    cell2csv([ '/home/ekaestne/PROJECTS/OUTPUT/mmps_compare/tract_FA/diffs' '/' mps_dev{1,trc_tfa_col_dev(iTF)} '_' mps_238{1,trc_tfa_col_238(iTF) } '.csv' ],avg_dff);
    
end

%% TRACT MD
ejk_chk_dir('/home/ekaestne/PROJECTS/OUTPUT/mmps_compare/tract_MD/figs')
ejk_chk_dir('/home/ekaestne/PROJECTS/OUTPUT/mmps_compare/tract_MD/diffs')

for iTF = 1:numel(trc_tmd_col_238)
    
    % Figure
    figure('Visible','off')
    scatter( cell2mat(mps_238(ovr_sbj_ind_238,trc_tmd_col_238(iTF))), cell2mat(mps_dev(ovr_sbj_ind_dev,trc_tmd_col_dev(iTF))), 'filled' )
    title( mmil_spec_char(mps_238{1,trc_tmd_col_238(iTF)},{'_' '-'},{' ' ' '}) )
    xlabel('MMPS 238')
    ylabel('MMPS DEV')
    print(gcf,[ '/home/ekaestne/PROJECTS/OUTPUT/mmps_compare/tract_MD/figs' '/' mps_dev{1,trc_tmd_col_dev(iTF)} '_' mps_238{1,trc_tmd_col_238(iTF) } '.png' ],'-dpng')
    close all
    
    % Diffs
    clear avg_dff
    for iS = 1:numel(ovr_sbj_ind_238)
        avg_dff{iS,1} = mps_238{ovr_sbj_ind_238(iS),2};
        avg_dff{iS,2} = mps_dev{ovr_sbj_ind_dev(iS),2};
        avg_dff{iS,3} = nanmean( cell2mat(mps_238(ovr_sbj_ind_238(iS),trc_tmd_col_238(iTF))) - ...
                                 cell2mat(mps_dev(ovr_sbj_ind_dev(iS),trc_tmd_col_dev(iTF))) );
    end
    cell2csv([ '/home/ekaestne/PROJECTS/OUTPUT/mmps_compare/tract_MD/diffs' '/' mps_dev{1,trc_tmd_col_dev(iTF)} '_' mps_238{1,trc_tmd_col_238(iTF) } '.csv' ],avg_dff);
    
end




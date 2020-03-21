clear; clc;

%%
prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'PhenotypeConnectome/';

grp_fle_hld = { 'Memory_Impairment_Connectome_SubjectList.csv' ... 'Language_Impairment_Connectome_SubjectList.csv' ...
                }; ...'Language_Impairment_Connectome_SubjectList_noOutliers.csv' };
grp_fle_nme = { 'all3T' ...
                }; ... 'subset' };
            
tmp_roi =  { { 'caudal-ITG'   'middle-ITG'         'rostral-ITG' ...
               'caudal-MTG'   'middle-MTG'         'rostral-MTG' ...
               'caudal-STG'   'middle-STG'         'rostral-STG' ...
               'temporalpole' 'transversetemporal' 'bankssts' ...
               'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform' ...
                'parahippocampal' 'entorhinal' } ...
           ...
           ...  { 'lh.caudal-ITG'      'lh.middle-ITG'         'lh.rostral-ITG' ...
           ...    'lh.caudal-MTG'      'lh.middle-MTG'         'lh.rostral-MTG' ...
           ...    'lh.caudal-STG'      'lh.middle-STG'         'lh.rostral-STG' ...
           ...    'lh.temporalpole'    'lh.transversetemporal' 'lh.bankssts' ...
           ...    'lh.caudal-fusiform' 'lh.middle-fusiform'    'lh.rostral-fusiform' ...
           ...    'lh.parahippocampal' 'lh.entorhinal' } ...
           ...
           ...  { 'caudal-MTG'   'middle-MTG'         'rostral-MTG' ...
           ...    'caudal-STG'   'middle-STG'         'rostral-STG' ...
           ...    'temporalpole' 'transversetemporal' 'bankssts' ...
           ...    'inferior-precentral' 'inferior-postcentral' ...
           ...    'parstriangularis'    'parsopercularis' ...
           ...    'supramarginal' 'insula' } ...
               ...
            ... { 'lh.caudal-MTG'   'lh.middle-MTG'         'lh.rostral-MTG' ...
            ...   'lh.caudal-STG'   'lh.middle-STG'         'lh.rostral-STG' ...
            ...   'lh.temporalpole' 'lh.transversetemporal' 'lh.bankssts' ...
            ...   'lh.inferior-precentral' 'lh.inferior-postcentral' ...
            ...   'lh.parstriangularis'    'lh.parsopercularis' ...
            ...   'lh.supramarginal' 'lh.insula' }  
             };
tmp_roi_nme = { 'temporal'        ...
                ... 'lhs_temporal'    ...
                ... 'perisylvian'     ...
               }; ... 'lhs_perisylvian' };
            
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iGR = 1:numel(grp_fle_nme)
    
    grp_fle = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' 'projects' '/' prj_nme '/' grp_fle_hld{iGR}]);
    grp_fle(cellfun(@isempty,grp_fle)) = {''};
    grp_fle(:,1) = strrep(grp_fle(:,1),'ducsf','d_ucsf');
    
    trn_tst = find(strcmpi(grp_fle(:,3),'UCSD'));
    tst_ind = find(strcmpi(grp_fle(:,3),'UCSF'));
    
    dta_lbl = double(strcmpi(grp_fle(:,2)      ,'memoryimpaired'));
    trn_lbl = double(strcmpi(grp_fle(trn_tst,2),'memoryimpaired'));
    tst_lbl = double(strcmpi(grp_fle(tst_ind,2),'memoryimpaired'));
    
    cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' grp_fle_nme{iGR} '_' 'totaldata_label' '.csv' ]  ,num2cell(dta_lbl))
    cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' grp_fle_nme{iGR} '_' 'traindata_label' '.csv' ] , num2cell(trn_lbl))
    cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' grp_fle_nme{iGR} '_' 'testdata_label' '.csv' ]  , num2cell(tst_lbl))
    
    for iTR = 1:numel(tmp_roi_nme)
        
        tmp_roi_hld = tmp_roi{iTR};
        
        % Load Data
        fcfg = [];
        fcfg.prj_dir = prj_dir;
        fcfg.sbj_nme = grp_fle(:,1);
        fcfg.dta_typ = 'dti_cnn';
        dti_cnn_dta = ejk_load_mcd_data(fcfg);
        
        top_tri = logical(triu(ones(size(dti_cnn_dta.dti_cnn_dta,2),size(dti_cnn_dta.dti_cnn_dta,3)),1));
        
        lbl_hld = dti_cnn_dta.cll_lbl(logical(top_tri));
        for iSB = 1:numel(dti_cnn_dta.sbj_nme)
            dta_hld = squeeze(dti_cnn_dta.dti_cnn_dta(iSB,:,:));
            dta_cnn(iSB,:) = dta_hld(top_tri)';
        end
        
        % ROI-to-All        
        lbl_spl_hld = regexp(lbl_hld,'==','split');
        hld_nme = [];
        for iT = 1:numel(tmp_roi_hld)
            hld_nme = [hld_nme find(~cellfun(@isempty,cellfun(@(x) string_find(x,tmp_roi_hld{iT}),lbl_spl_hld,'uni',0)))'];
        end
        hld_nme = unique(hld_nme);
        
        dta_cnn_tmp = dta_cnn(:,hld_nme);
        lbl_hld_tmp = lbl_hld(hld_nme);
        
        dta_trn_cnn_tmp = dta_cnn(trn_tst,:);
        dta_tst_cnn_tmp = dta_cnn(tst_ind,:);
        
        cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' grp_fle_nme{iGR} '_' 'allROI' '_' 'toall' '_' 'totaldata' '.csv' ] ,num2cell(dta_cnn_tmp))
        cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' grp_fle_nme{iGR} '_' 'allROI' '_' 'toall' '_' 'traindata' '.csv' ] ,num2cell(dta_trn_cnn_tmp))
        cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' grp_fle_nme{iGR} '_' 'allROI' '_' 'toall' '_' 'testdata' '.csv' ]  ,num2cell(dta_tst_cnn_tmp))
        
        % ROI-to-ROI
        pos_hld = zeros(numel(lbl_spl_hld),2);
        for iPS = 1:numel(lbl_spl_hld)
            pos_hld(iPS,1) = any(~cellfun(@isempty,cellfun(@(x) string_find( lbl_spl_hld{iPS}(1) , x ) , tmp_roi_hld , 'uni', 0)));
            pos_hld(iPS,2) =  any(~cellfun(@isempty,cellfun(@(x) string_find( lbl_spl_hld{iPS}(2) , x ) , tmp_roi_hld , 'uni', 0)));
        end
        hld_nme = find(sum(pos_hld,2)==2);
        
        dta_cnn_tmp = dta_cnn(:,hld_nme);
        lbl_hld_tmp = lbl_hld(hld_nme);
        
        dta_trn_cnn_tmp = dta_cnn_tmp(trn_tst,:);
        dta_tst_cnn_tmp = dta_cnn_tmp(tst_ind,:);
        
        cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' grp_fle_nme{iGR} '_' tmp_roi_nme{iTR} '_' 'only' '_' 'totaldata' '.csv' ] ,num2cell(dta_cnn_tmp))
        cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' grp_fle_nme{iGR} '_' tmp_roi_nme{iTR} '_' 'only' '_' 'traindata' '.csv' ] ,num2cell(dta_trn_cnn_tmp))
        cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' grp_fle_nme{iGR} '_' tmp_roi_nme{iTR} '_' 'only' '_' 'testdata' '.csv' ]  ,num2cell(dta_tst_cnn_tmp))       
        
    end
end

%%


























%% Original Subject List - Temporal Subnet
% Loading 1a) ALL EPD Subjects
grp_fle = 'Language_Impairment_Connectome_SubjectList.csv';

grp_fle = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' 'tmp' '/' prj_nme '/' grp_fle]);
    grp_fle(cellfun(@isempty,grp_fle)) = {''};
    
% Load Connectome Data
fcfg = [];
fcfg.prj_dir = prj_dir;
fcfg.sbj_nme = grp_fle(:,1);
fcfg.dta_typ = 'dti_cnn';
dti_cnn_dta = ejk_load_mcd_data(fcfg);

top_tri = logical(triu(ones(size(dti_cnn_dta.dti_cnn_dta,2),size(dti_cnn_dta.dti_cnn_dta,3)),1));

lbl_hld = dti_cnn_dta.cll_lbl(logical(top_tri));
for iSB = 1:numel(dti_cnn_dta.sbj_nme)
    dta_hld = squeeze(dti_cnn_dta.dti_cnn_dta(iSB,:,:));
    dta_cnn(iSB,:) = dta_hld(top_tri)';
end

% 2) Get Temporal Subnet
tmp_roi = { 'caudal-ITG'   'middle-ITG'         'rostral-ITG' ...
            'caudal-MTG'   'middle-MTG'         'rostral-MTG' ...
            'caudal-STG'   'middle-STG'         'rostral-STG' ...
            'temporalpole' 'transversetemporal' 'bankssts' ...
            'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform' ...
            'parahippocampal' 'entorhinal' };

lbl_spl_hld = regexp(lbl_hld,'==','split');
hld_nme = [];
for iTR = 1:numel(tmp_roi)
    hld_nme = [hld_nme find(~cellfun(@isempty,cellfun(@(x) string_find(x,tmp_roi{iTR}),lbl_spl_hld,'uni',0)))'];
end
hld_nme = unique(hld_nme);

dta_cnn_tmp = dta_cnn(:,hld_nme);
lbl_hld_tmp = lbl_hld(hld_nme);

% Split into 3) UCSD TRAIN/TEST
trn_tst = find(strcmpi(grp_fle(:,3),'UCSD'));
tst_ind = find(strcmpi(grp_fle(:,3),'UCSF'));

dta_lbl = double(strcmpi(grp_fle(:,2),'languageimpaired'));
trn_lbl = double(strcmpi(grp_fle(trn_tst,2),'languageimpaired'));
tst_lbl = double(strcmpi(grp_fle(tst_ind,2),'languageimpaired'));

dta_trn_cnn_tmp = dta_cnn_tmp(trn_tst,:);
dta_tst_cnn_tmp = dta_cnn_tmp(tst_ind,:);

% Save
cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'connectionname' '.csv' ] ,lbl_hld_tmp)
    
cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'totaldata' '.csv' ]            ,num2cell(dta_cnn_tmp))
    cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'totaldata_label' '.csv' ]  ,num2cell(dta_lbl))

cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'traindata' '.csv' ]           ,num2cell(dta_trn_cnn_tmp))
    cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'traindata_label' '.csv' ] ,num2cell(trn_lbl))
cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'testdata' '.csv' ]            ,num2cell(dta_tst_cnn_tmp))
    cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'testdata_label' '.csv' ]  ,num2cell(tst_lbl))

%% Cleand Subject List - Temporal Subnet
clear dta_cnn

% Loading 1b) Cleaned EPD Subjects
grp_fle = 'Language_Impairment_Connectome_SubjectList_noOutliers.csv';

grp_fle = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' 'tmp' '/' prj_nme '/' grp_fle]);
    grp_fle(cellfun(@isempty,grp_fle)) = {''};
    
% Load Connectome Data
fcfg = [];
fcfg.prj_dir = prj_dir;
fcfg.sbj_nme = grp_fle(:,1);
fcfg.dta_typ = 'dti_cnn';
dti_cnn_dta = ejk_load_mcd_data(fcfg);

top_tri = logical(triu(ones(size(dti_cnn_dta.dti_cnn_dta,2),size(dti_cnn_dta.dti_cnn_dta,3)),1));

lbl_hld = dti_cnn_dta.cll_lbl(logical(top_tri));
for iSB = 1:numel(dti_cnn_dta.sbj_nme)
    dta_hld = squeeze(dti_cnn_dta.dti_cnn_dta(iSB,:,:));
    dta_cnn(iSB,:) = dta_hld(top_tri)';
end

% 2) Get Temporal Subnet
tmp_roi = { 'caudal-ITG'   'middle-ITG'         'rostral-ITG' ...
            'caudal-MTG'   'middle-MTG'         'rostral-MTG' ...
            'caudal-STG'   'middle-STG'         'rostral-STG' ...
            'temporalpole' 'transversetemporal' 'bankssts' ...
            'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform' ...
            'parahippocampal' 'entorhinal' };

lbl_spl_hld = regexp(lbl_hld,'==','split');
hld_nme = [];
for iTR = 1:numel(tmp_roi)
    hld_nme = [hld_nme find(~cellfun(@isempty,cellfun(@(x) string_find(x,tmp_roi{iTR}),lbl_spl_hld,'uni',0)))'];
end
hld_nme = unique(hld_nme);

dta_cnn_tmp = dta_cnn(:,hld_nme);
lbl_hld_tmp = lbl_hld(hld_nme);

% Split into 3) UCSD TRAIN/TEST
trn_tst = find(strcmpi(grp_fle(:,3),'UCSD'));
tst_ind = find(strcmpi(grp_fle(:,3),'UCSF'));
con_ind = find(strcmpi(grp_fle(:,3),'HC'));

dta_lbl = double(strcmpi(grp_fle(:,2),'languageimpaired'));
trn_lbl = double(strcmpi(grp_fle(trn_tst,2),'languageimpaired'));

dta_trn_cnn_tmp = dta_cnn_tmp(trn_tst,:);
dta_tst_cnn_tmp = dta_cnn_tmp(tst_ind,:);
dta_con_cnn_tmp = dta_cnn_tmp(con_ind,:);

% Save
cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'NoOutliers' '_' 'connectionname' '.csv' ] ,lbl_hld_tmp)

cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'NoOutliers' '_' 'totaldata' '.csv' ]            ,num2cell(dta_cnn_tmp))
    cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'NoOutliers' '_' 'totaldata_label' '.csv' ]  ,num2cell(dta_lbl))
    
cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'NoOutliers' '_' 'traindata' '.csv' ]           ,num2cell(dta_trn_cnn_tmp))
    cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'NoOutliers' '_' 'traindata_label' '.csv' ] ,num2cell(trn_lbl))
cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'NoOutliers' '_' 'testdata' '.csv' ]            ,num2cell(dta_tst_cnn_tmp))
    cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'NoOutliers' '_' 'testdata_label' '.csv' ]  ,num2cell(tst_lbl))

%% ALL Subject List - Left Temporal Subnet
% Loading 1b) Cleaned EPD Subjects
grp_fle = 'Language_Impairment_Connectome_SubjectList.csv';

grp_fle = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' 'tmp' '/' prj_nme '/' grp_fle]);
    grp_fle(cellfun(@isempty,grp_fle)) = {''};
    
% Load Connectome Data
fcfg = [];
fcfg.prj_dir = prj_dir;
fcfg.sbj_nme = grp_fle(:,1);
fcfg.dta_typ = 'dti_cnn';
dti_cnn_dta = ejk_load_mcd_data(fcfg);

top_tri = logical(triu(ones(size(dti_cnn_dta.dti_cnn_dta,2),size(dti_cnn_dta.dti_cnn_dta,3)),1));

lbl_hld = dti_cnn_dta.cll_lbl(logical(top_tri));
for iSB = 1:numel(dti_cnn_dta.sbj_nme)
    dta_hld = squeeze(dti_cnn_dta.dti_cnn_dta(iSB,:,:));
    dta_cnn(iSB,:) = dta_hld(top_tri)';
end

% 2) Get Temporal Subnet
tmp_roi = { 'lh.caudal-ITG'      'lh.middle-ITG'         'lh.rostral-ITG' ...
            'lh.caudal-MTG'      'lh.middle-MTG'         'lh.rostral-MTG' ...
            'lh.caudal-STG'      'lh.middle-STG'         'lh.rostral-STG' ...
            'lh.temporalpole'    'lh.transversetemporal' 'lh.bankssts' ...
            'lh.caudal-fusiform' 'lh.middle-fusiform'    'lh.rostral-fusiform' ...
            'lh.parahippocampal' 'lh.entorhinal' };

lbl_spl_hld = regexp(lbl_hld,'==','split');
hld_nme = [];
for iTR = 1:numel(tmp_roi)
    hld_nme = [hld_nme find(~cellfun(@isempty,cellfun(@(x) string_find(x,tmp_roi{iTR}),lbl_spl_hld,'uni',0)))'];
end
hld_nme = unique(hld_nme);

dta_cnn_tmp = dta_cnn(:,hld_nme);
lbl_hld_tmp = lbl_hld(hld_nme);

% Split into 3) UCSD TRAIN/TEST
trn_tst = find(strcmpi(grp_fle(:,3),'UCSD'));
tst_ind = find(strcmpi(grp_fle(:,3),'UCSF'));
con_ind = find(strcmpi(grp_fle(:,3),'HC'));

dta_lbl = double(strcmpi(grp_fle(:,2),'languageimpaired'));
trn_lbl = double(strcmpi(grp_fle(trn_tst,2),'languageimpaired'));

dta_trn_cnn_tmp = dta_cnn_tmp(trn_tst,:);
dta_tst_cnn_tmp = dta_cnn_tmp(tst_ind,:);
dta_con_cnn_tmp = dta_cnn_tmp(con_ind,:);

% Save
cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'connectionname'  '_leftonly.csv' ] ,lbl_hld_tmp)

cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'totaldata'  '_leftonly.csv' ]            ,num2cell(dta_cnn_tmp))
    cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'totaldata_label'  '_leftonly.csv' ]  ,num2cell(dta_lbl))

cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'condata'  '_leftonly.csv' ] , num2cell(dta_con_cnn_tmp) )
    
cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'traindata'  '_leftonly.csv' ]           ,num2cell(dta_trn_cnn_tmp))
    cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'traindata_label'  '_leftonly.csv' ] ,num2cell(trn_lbl))
cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'testdata'  '_leftonly.csv' ]            ,num2cell(dta_tst_cnn_tmp))
    cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'testdata_label'  '_leftonly.csv' ]  ,num2cell(tst_lbl))

    %% ALL Subject List - Left Temporal Subnet
% Loading 1b) Cleaned EPD Subjects
grp_fle = 'Language_Impairment_Connectome_SubjectList.csv';

grp_fle = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' 'tmp' '/' prj_nme '/' grp_fle]);
    grp_fle(cellfun(@isempty,grp_fle)) = {''};
    
% Load Connectome Data
fcfg = [];
fcfg.prj_dir = prj_dir;
fcfg.sbj_nme = grp_fle(:,1);
fcfg.dta_typ = 'dti_cnn';
dti_cnn_dta = ejk_load_mcd_data(fcfg);

top_tri = logical(triu(ones(size(dti_cnn_dta.dti_cnn_dta,2),size(dti_cnn_dta.dti_cnn_dta,3)),1));

lbl_hld = dti_cnn_dta.cll_lbl(logical(top_tri));
for iSB = 1:numel(dti_cnn_dta.sbj_nme)
    dta_hld = squeeze(dti_cnn_dta.dti_cnn_dta(iSB,:,:));
    dta_cnn(iSB,:) = dta_hld(top_tri)';
end

% 2) Get Temporal Subnet
tmp_roi = { 'lh.caudal-ITG'      'lh.middle-ITG'         'lh.rostral-ITG' ...
            'lh.caudal-MTG'      'lh.middle-MTG'         'lh.rostral-MTG' ...
            'lh.caudal-STG'      'lh.middle-STG'         'lh.rostral-STG' ...
            'lh.temporalpole'    'lh.transversetemporal' 'lh.bankssts' ...
            'lh.caudal-fusiform' 'lh.middle-fusiform'    'lh.rostral-fusiform' ...
            'lh.parahippocampal' 'lh.entorhinal' };

lbl_spl_hld = regexp(lbl_hld,'==','split');
hld_nme = [];
for iTR = 1:numel(tmp_roi)
    hld_nme = [hld_nme find(~cellfun(@isempty,cellfun(@(x) string_find(x,tmp_roi{iTR}),lbl_spl_hld,'uni',0)))'];
end
hld_nme = unique(hld_nme);

dta_cnn_tmp = dta_cnn(:,hld_nme);
lbl_hld_tmp = lbl_hld(hld_nme);

% Split into 3) UCSD TRAIN/TEST
trn_tst = find(strcmpi(grp_fle(:,3),'UCSD'));
tst_ind = find(strcmpi(grp_fle(:,3),'UCSF'));
con_ind = find(strcmpi(grp_fle(:,3),'HC'));

dta_lbl = double(strcmpi(grp_fle(:,2),'languageimpaired'));
trn_lbl = double(strcmpi(grp_fle(trn_tst,2),'languageimpaired'));

dta_trn_cnn_tmp = dta_cnn_tmp(trn_tst,:);
dta_tst_cnn_tmp = dta_cnn_tmp(tst_ind,:);
dta_con_cnn_tmp = dta_cnn_tmp(con_ind,:);

% Save
cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'connectionname'  '_leftonly.csv' ] ,lbl_hld_tmp)

cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'totaldata'  '_leftonly.csv' ]            ,num2cell(dta_cnn_tmp))
    cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'totaldata_label'  '_leftonly.csv' ]  ,num2cell(dta_lbl))

cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'condata'  '_leftonly.csv' ] , num2cell(dta_con_cnn_tmp) )
    
cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'traindata'  '_leftonly.csv' ]           ,num2cell(dta_trn_cnn_tmp))
    cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'traindata_label'  '_leftonly.csv' ] ,num2cell(trn_lbl))
cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'testdata'  '_leftonly.csv' ]            ,num2cell(dta_tst_cnn_tmp))
    cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'testdata_label'  '_leftonly.csv' ]  ,num2cell(tst_lbl))

%% ALL Subject List - Perisylvian
% Loading 1b) Cleaned EPD Subjects
grp_fle = 'Language_Impairment_Connectome_SubjectList.csv';

grp_fle = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' 'tmp' '/' prj_nme '/' grp_fle]);
    grp_fle(cellfun(@isempty,grp_fle)) = {''};
    
% Load Connectome Data
fcfg = [];
fcfg.prj_dir = prj_dir;
fcfg.sbj_nme = grp_fle(:,1);
fcfg.dta_typ = 'dti_cnn';
dti_cnn_dta = ejk_load_mcd_data(fcfg);

top_tri = logical(triu(ones(size(dti_cnn_dta.dti_cnn_dta,2),size(dti_cnn_dta.dti_cnn_dta,3)),1));

lbl_hld = dti_cnn_dta.cll_lbl(logical(top_tri));
for iSB = 1:numel(dti_cnn_dta.sbj_nme)
    dta_hld = squeeze(dti_cnn_dta.dti_cnn_dta(iSB,:,:));
    dta_cnn(iSB,:) = dta_hld(top_tri)';
end

% 2) Get Temporal Subnet
tmp_roi = { 'caudal-MTG'   'middle-MTG'         'rostral-MTG' ...
            'caudal-STG'   'middle-STG'         'rostral-STG' ...
            'temporalpole' 'transversetemporal' 'bankssts' ...
            'inferior-precentral' 'inferior-postcentral' ...
            'parstriangularis'    'parsopercularis' ...
            'supramarginal' 'insula' };

lbl_spl_hld = regexp(lbl_hld,'==','split');
hld_nme = [];
for iTR = 1:numel(tmp_roi)
    hld_nme = [hld_nme find(~cellfun(@isempty,cellfun(@(x) string_find(x,tmp_roi{iTR}),lbl_spl_hld,'uni',0)))'];
end
hld_nme = unique(hld_nme);

dta_cnn_tmp = dta_cnn(:,hld_nme);
lbl_hld_tmp = lbl_hld(hld_nme);

% Split into 3) UCSD TRAIN/TEST
trn_tst = find(strcmpi(grp_fle(:,3),'UCSD'));
tst_ind = find(strcmpi(grp_fle(:,3),'UCSF'));
con_ind = find(strcmpi(grp_fle(:,3),'HC'));

dta_lbl = double(strcmpi(grp_fle(:,2),'languageimpaired'));
trn_lbl = double(strcmpi(grp_fle(trn_tst,2),'languageimpaired'));
tst_lbl = double(strcmpi(grp_fle(tst_ind,2),'languageimpaired'));

dta_trn_cnn_tmp = dta_cnn_tmp(trn_tst,:);
dta_tst_cnn_tmp = dta_cnn_tmp(tst_ind,:);

% Save
cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'connectionname'  '_peri.csv' ] ,lbl_hld_tmp)

cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'totaldata'  '_peri.csv' ]            ,num2cell(dta_cnn_tmp))
    cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'totaldata_label'  '_peri.csv' ]  ,num2cell(dta_lbl))
   
cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'traindata'  '_peri.csv' ]           ,num2cell(dta_trn_cnn_tmp))
    cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'traindata_label'  '_peri.csv' ] ,num2cell(trn_lbl))
cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'testdata'  '_peri.csv' ]            ,num2cell(dta_tst_cnn_tmp))
    cell2csv(['/home/ekaestne/PROJECTS/DATA/memory_connectome' '/' 'Data' '_' 'testdata_label'  '_peri.csv' ]  ,num2cell(tst_lbl))

    
    
    
    
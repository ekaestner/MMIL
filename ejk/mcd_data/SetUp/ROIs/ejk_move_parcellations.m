%
%
%
%
%

function ejk_move_parcellations

hme_dir = '/home/ekaestne/PROJECTS';

roi_nme = { 'aparc.' ...
            'aparc.a2009s.' };

%%
sbj_nme = mmil_readtext( [ hme_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_subject_name.csv' ] );
    sbj_nme{end+1} = 'fsaverage';
fsr_nme = mmil_readtext( [ hme_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_freesurfer_recons.csv' ] );
    fsr_nme{end+1,1} = 'fsaverage';
    fsr_nme{end,2}   = '/home/mmilmcdRSI/data/fsurf';
    fsr_nme{end,3}   = 'fsaverage';
    
%% Move over Desikan & Destrieaux
for iSB = 1:size( sbj_nme , 1 )
    
    if ~isempty(fsr_nme{iSB,3})
     
        if ~strcmpi( fsr_nme{iSB,1} , 'fsaverage' )
            sbj_dir_lst = dir(sprintf('%s/FSURF_*',fsr_nme{iSB,2}));
            sbj_dir = regexp({sbj_dir_lst.name},['FSURF_' fsr_nme{iSB,3} '.+_1$'],'match'); sbj_dir = [sbj_dir{:}];
        else
            sbj_dir_lst = dir(sprintf('%s/*',fsr_nme{iSB,2}));
            sbj_dir = regexp({sbj_dir_lst.name},['' fsr_nme{iSB,3}],'match'); sbj_dir = [sbj_dir{:}];
        end
    
    ejk_chk_dir( [ hme_dir '/' 'DATA' '/' fsr_nme{iSB,1} '/' ] );
    ejk_chk_dir( [ hme_dir '/' 'DATA' '/' fsr_nme{iSB,1} '/' 'ROIs' '/' ] );
    
    for iRN = 1:numel(roi_nme)
        if ~exist( [ hme_dir '/' 'DATA' '/' fsr_nme{iSB,1} '/' 'ROIs' '/' 'lh.' roi_nme{iRN} 'annot' ] , 'file' ) && ...
            ~isempty(sbj_dir)
            
            copyfile( [ fsr_nme{iSB,2} '/' sbj_dir{:} '/' 'label' '/' 'lh.' roi_nme{iRN} 'annot' ] , ...
                      [ hme_dir '/' 'DATA' '/' fsr_nme{iSB,1} '/' 'ROIs' '/' 'lh.' roi_nme{iRN} 'annot' ] );
            copyfile( [ fsr_nme{iSB,2} '/' sbj_dir{:} '/' 'label' '/' 'rh.' roi_nme{iRN} 'annot' ] , ...
                      [ hme_dir '/' 'DATA' '/' fsr_nme{iSB,1} '/' 'ROIs' '/' 'rh.' roi_nme{iRN} 'annot' ] );
        end
    end
    end
end

%% Create Modified Desikan
for iSB = 1:size( sbj_nme , 1 )
    
    fcfg = [];
    
    fcfg.prj_dir = hme_dir;
    fcfg.sbj_nme = sbj_nme{iSB};
   
    fcfg.fsr_dir = fsr_nme{iSB,2};
    fcfg.fsr_nme = fsr_nme{iSB,3};

    fcfg.ovr_wrt = 1;

    mmil_split_parcellation(fcfg)
        
end

%% Create Additional ROIs


end




















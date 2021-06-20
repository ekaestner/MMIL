clear; clc;

emy_cph = mmil_readtext([ '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Emory_Slah_Explore' '/' 'subject_cipher.csv' ]);

ntl_dta = mmil_readtext([ '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Emory_Slah_Explore' '/' 'NatalieData.csv' ]);

prc_dta = mmil_readtext([ '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Emory_Slah_Explore/EmorySlahExplore_v2/4th' '/' 'Emory_Scan_Classifications.csv' ]);
    prc_dta = prc_dta(2:end,:);

%%
ejk_dta_out = cell( size(prc_dta,1), 5 );



for iS = 1:size(prc_dta,1)
  
    ejk_dta_out{iS,2} = prc_dta{iS,1};
    cph_ind = find(strcmpi( emy_cph(:,4), prc_dta{iS,1} ));
    
    if ~isempty(cph_ind)
        sbj_nme = emy_cph{cph_ind,1};
        emy_nme = emy_cph{cph_ind,2};
        
        ntl_ind = find(strcmpi( ntl_dta(:,2), emy_nme));
        
        ejk_dta_out{iS,1} = sbj_nme;
        ejk_dta_out{iS,2} = prc_dta{iS,1};
        ejk_dta_out{iS,3} = emy_nme;
        
        ejk_dta_out{iS,4} = ntl_dta{ ntl_ind, 2};
        ejk_dta_out{iS,5} = ntl_dta{ ntl_ind, 5};
        ejk_dta_out{iS,6} = ntl_dta{ ntl_ind, 4};
    end
    
end

cell2csv( ['/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/slh_atl_mem' '/' 'sbj_nme.csv'] , ejk_dta_out);

%%
row_nme = { 'Date of Scan' 'Date of ' ... 
            'LM I' 'LM II' 'VPA I' 'VPA II' 'BVMT Immediate' 'BVMT Delayed' ...
            'Age' 'Education' 'Sex (F/M)' 'Handedness (R/L)' 'Race (White/>1/Black/Asian/NH/Unknown' 'Ethnicity (Non-Hispanic/Hispanic/Unknown' ...
            'Premorbid estimated IQ (WTAR)' 'Age of Seizure Onset' 'Duration of Epilepsy' 'Number of ASMs' 'MTS (Y/N)' ...
            'Resection Side (Dominant/Non-dominant' 'Seizure frequency' 'Lifetime GTC frequency (0-1/2-9/10-39/>40' ...
            'Engel Outcome (I/II+)' };

csv_out = [ {'SubjID' 'FreesurferRecon' 'Surgery' 'Side Of Seizure Onset'} row_nme ; ejk_dta_out(:,[1 2 5 6])  cell(size(ejk_dta_out,1),size(row_nme,2)) ];
cell2csv( ['/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/slh_atl_mem' '/' 'Emory_Data_For_UCSD.csv'] , csv_out);
























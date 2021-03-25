clear; clc;

prj_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/CoAuthor/Natalie/';
out_dir = '/home/ekaestne/PROJECTS/DATA/Natalie_BNT_Project';

%%
sbj_nme = mmil_readtext([prj_dir '/' 'ForNatalie_201027.csv']);
sbj_nme_col = sbj_nme(1,:);
sbj_nme     = sbj_nme(2:end, :);

dta_loc     = mmil_readtext([prj_dir '/' 'dti_search_noapostrophe.csv']);
sbj_dbl_mri =  mmil_readtext([out_dir '/' 'Multiple_FSURFs_mri_fixed.csv']);
sbj_dbl_dti =  mmil_readtext([out_dir '/' 'Multiple_FSURFs_dti_fixed.csv']);

dti_sbj_ind = strcmpi( sbj_nme( :, find(strcmpi(sbj_nme_col, 'DTI available?')) ), 'yes');
t3_sbj_ind  = cell2mat(sbj_nme( :, find(strcmpi(sbj_nme_col, 'Field Strength')) )) == 3;
bnt_sbj_ind = strcmpi( sbj_nme( :, find(strcmpi(sbj_nme_col, 'BNT available?')) ), 'yes');
% Exclude B0 scans

sbj_nme_use = sbj_nme( dti_sbj_ind & t3_sbj_ind & bnt_sbj_ind, 1);

%% DTI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
for iS = 1:size(sbj_nme_use,1)
    row_ind = string_find( dta_loc(:,1),sbj_nme_use{iS} );
    loc_hld{iS} = cat( 1, dta_loc(row_ind, [1 2 4 7 ] ) );
end

loc_ind = cell( size(sbj_nme_use,1), 5 );
for iS = 1:size(sbj_nme_use,1)
    
    san_die_loc = '30-dir DTI Pepolar #1';
    san_fra_loc = '(B=1000 DIR=30)';
    
    loc_ind{iS,1} = sbj_nme_use{iS,1};
    loc_ind{iS,2} = string_find( loc_hld{iS}(:,4), san_die_loc );
    loc_ind{iS,3} = string_find( loc_hld{iS}(:,4), san_fra_loc );
    loc_ind{iS,4} = sum( [numel(loc_ind{iS,2}) numel(loc_ind{iS,3})] ) > 0;
    loc_ind{iS,5} = [ loc_ind{iS,2} loc_ind{iS,3} ];
    
end

%%
out_tbl = cell( size(sbj_nme_use,1), 15 );
dub_tbl = cell( size(sbj_nme_use,1), 2 );

grp_col_ind     = find( strcmpi( sbj_nme_col, 'Group'));
age_ons_col_ind = find( strcmpi( sbj_nme_col, 'Age of Epilepsy Onset (yrs)'));
sde_ons_col_ind = find( strcmpi( sbj_nme_col, 'Side of Epileptic Focus'));
mst_col_ind     = find( strcmpi( sbj_nme_col, 'Does the patient have Mesial Temporal Pathology (e.g. MTS)'));

age_col_ind     = find( strcmpi( sbj_nme_col, 'Age at Testing'));
hnd_col_ind     = find( strcmpi( sbj_nme_col, 'Handedness'));
edu_col_ind     = find( strcmpi( sbj_nme_col, 'Years of Education'));
lng_grp_col_ind = find( strcmpi( sbj_nme_col, 'Language Group'));
gdr_col_ind     = find( strcmpi( sbj_nme_col, 'Gender'));

bnt_scr_col_ind = find( strcmpi( sbj_nme_col, 'Boston Naming Test  Total Correct'));
bnt_tsc_col_ind = find( strcmpi( sbj_nme_col, 'Boston Naming Test Tscore'));

fld_str_col_ind = find( strcmpi( sbj_nme_col, 'Field Strength'));

out_csv_nme = { 'SubjectID' ...
    sbj_nme_col{grp_col_ind} ...
    sbj_nme_col{age_ons_col_ind} ...
    sbj_nme_col{sde_ons_col_ind} ...
    sbj_nme_col{mst_col_ind} ...
    sbj_nme_col{age_col_ind} ...
    sbj_nme_col{hnd_col_ind} ...
    sbj_nme_col{edu_col_ind} ...
    sbj_nme_col{lng_grp_col_ind} ...
    sbj_nme_col{gdr_col_ind} ...
    sbj_nme_col{bnt_scr_col_ind} ...
    sbj_nme_col{bnt_tsc_col_ind} ...
    sbj_nme_col{fld_str_col_ind} };

for iS = 1:size(sbj_nme_use,1)
    
    sbj_nme_ind = find( strcmpi( sbj_nme(:,1), sbj_nme_use{iS,1}) );
    
    out_tbl{iS,1} = sbj_nme_use{iS, 1};
    
    out_tbl{iS,2} = sbj_nme{sbj_nme_ind, grp_col_ind };
    out_tbl{iS,3} = sbj_nme{sbj_nme_ind, age_ons_col_ind };
    out_tbl{iS,4} = sbj_nme{sbj_nme_ind, sde_ons_col_ind };
    out_tbl{iS,5} = sbj_nme{sbj_nme_ind, mst_col_ind };
    
    out_tbl{iS,6}  = sbj_nme{sbj_nme_ind, age_col_ind };
    out_tbl{iS,7}  = sbj_nme{sbj_nme_ind, hnd_col_ind };
    out_tbl{iS,8}  = sbj_nme{sbj_nme_ind, edu_col_ind };
    out_tbl{iS,9}  = sbj_nme{sbj_nme_ind, lng_grp_col_ind };
    out_tbl{iS,10} = sbj_nme{sbj_nme_ind, gdr_col_ind };
    
    out_tbl{iS,11} = sbj_nme{sbj_nme_ind, bnt_scr_col_ind };
    out_tbl{iS,12} = sbj_nme{sbj_nme_ind, bnt_tsc_col_ind };
    
    out_tbl{iS,13} = sbj_nme{sbj_nme_ind, fld_str_col_ind };
    
    if numel(loc_ind{iS,5}) == 0
        out_tbl{iS,14}  = 'Empty';
        out_tbl{iS,15} = 'N/A';
        out_tbl{iS,16} = 'N/A';
    elseif numel(loc_ind{iS,5}) > 1
        
        dbl_ind = find( strcmpi( sbj_dbl_dti(:,1), sbj_nme_use{iS,1} ));
        
        if ~isempty(dbl_ind)
            chs_hld = string_find( loc_hld{iS}(:,1), sbj_dbl_dti{ dbl_ind, 2} );
            chs_loc = intersect( loc_ind{iS,5}, chs_hld);
            
            end_pnt = strfind( loc_hld{iS}{ chs_loc , 1 }, '/')-1;
            out_tbl{iS,14}  = loc_hld{iS}{ chs_loc , 1 }(1:end_pnt);
            out_tbl{iS,15} = loc_hld{iS}{ chs_loc , 2 };
            out_tbl{iS,16} = loc_hld{iS}{ chs_loc , 4 };
            
        else
            out_tbl{iS,14}  = 'Multiple';
            out_tbl{iS,15} = 'N/A';
            out_tbl{iS,16} = 'N/A';
            
            dub_tbl{iS,1} = sbj_nme_use{iS,1};
            dub_tbl{iS,2} = strcat( loc_hld{iS}( loc_ind{iS,5} , 1), '; ');
            dub_tbl{iS,2} = [dub_tbl{iS,2}{:}];
        end
        
    else
        end_pnt        = strfind( loc_hld{iS}{ loc_ind{iS,5} , 1 }, '/')-1;
        out_tbl{iS,14}  = loc_hld{iS}{ loc_ind{iS,5} , 1 }(1:end_pnt);
        out_tbl{iS,15} = loc_hld{iS}{ loc_ind{iS,5} , 2 };
        out_tbl{iS,16} = loc_hld{iS}{ loc_ind{iS,5} , 4 };
    end
    
end

dub_tbl( cellfun( @isempty, dub_tbl(:,1)), :) = [];
cell2csv([ prj_dir '/' 'Multiple_FSURFs_dti.csv'],dub_tbl);

mss_tbl = out_tbl;
mss_tbl( ~cellfun( @isempty, mss_tbl(:,2)), :) = [];

hav_tbl = out_tbl;
hav_tbl( cellfun( @isempty, hav_tbl(:,2)), :)   = [];
hav_tbl( strcmpi( hav_tbl(:,14),'Multiple'), :) = [];

prb_tbl = hav_tbl;
prb_tbl( ~strcmpi(prb_tbl(:,14), 'N/A'), :) = [];

hav_tbl( strcmpi(hav_tbl(:,14), 'N/A'), :) = [];
cell2csv([ prj_dir '/' 'FSURFs_dti.csv'],hav_tbl);

[ size(sbj_nme_use,1) size(hav_tbl,1) size(mss_tbl,1) size(prb_tbl,1) size(dub_tbl,1) ]

%% Creat csv
col_nme = out_csv_nme;
cell2csv( [ out_dir '/' 'Subject_data.csv'], [ col_nme ; hav_tbl(:, 1:13) ]);

%% Create commands
frs_dir = '/space/syn09/1/data/MMILDB/MCD_RSI/raw/';
sec_dir = out_dir;

prn_cmd_dti = '';
for iS = 1:size(hav_tbl,1)
    
    prn_cmd_dti = sprintf('%s\n\n########################\n\nmkdir %s/%s/\n', prn_cmd_dti, sec_dir, hav_tbl{iS,1});
    prn_cmd_dti = sprintf('%smkdir %s/%s/dti/\n', prn_cmd_dti, sec_dir, hav_tbl{iS,1});
    prn_cmd_dti = sprintf('%scp -rL %s/%s/%s/*.dcm %s/%s/dti/\n', prn_cmd_dti, frs_dir, hav_tbl{iS,14}, hav_tbl{iS,15}, sec_dir, hav_tbl{iS,1});
    prn_cmd_dti = sprintf('%sls -l %s/%s/dti/ | wc -l\n\n\n\n', prn_cmd_dti, sec_dir, hav_tbl{iS,1});
    
end

cell2csv( [ out_dir '/' 'dti_copy_commands.csv'], {prn_cmd_dti});

%% MRI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

prj_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/CoAuthor/Natalie/';
out_dir = '/home/ekaestne/PROJECTS/DATA/Natalie_BNT_Project';

%%
sbj_nme = mmil_readtext([prj_dir '/' 'ForNatalie_201027.csv']);
sbj_nme_col = sbj_nme(1,:);
sbj_nme     = sbj_nme(2:end, :);

sbj_dbl_mri =  mmil_readtext([out_dir '/' 'Multiple_FSURFs_mri_fixed.csv']);
sbj_dbl_dti =  mmil_readtext([out_dir '/' 'Multiple_FSURFs_dti_fixed.csv']);

dti_sbj_ind = strcmpi( sbj_nme( :, find(strcmpi(sbj_nme_col, 'DTI available?')) ), 'yes');
t3_sbj_ind  = cell2mat(sbj_nme( :, find(strcmpi(sbj_nme_col, 'Field Strength')) )) == 3;
bnt_sbj_ind = strcmpi( sbj_nme( :, find(strcmpi(sbj_nme_col, 'BNT available?')) ), 'yes');
% Exclude B0 scans

sbj_nme_use = sbj_nme( dti_sbj_ind & t3_sbj_ind & bnt_sbj_ind, 1);

%%
dta_loc = mmil_readtext([prj_dir '/' 'mri_search_noapostrophe.csv']);
for iS = 1:size(sbj_nme_use,1)
    row_ind = string_find( dta_loc(:,1),sbj_nme_use{iS} );
    loc_hld{iS} = cat( 1, dta_loc(row_ind, [1 2 4 7 ] ) );
end

loc_ind = cell( size(sbj_nme_use,1), 5 );
for iS = 1:size(sbj_nme_use,1)
    
    san_die_loc = 'FSPGR_SAG_';
    san_fra_loc = 'IRSPGR_';
    
    loc_ind{iS,1} = sbj_nme_use{iS,1};
    loc_ind{iS,2} = string_find( loc_hld{iS}(:,4), san_die_loc );
    loc_ind{iS,3} = string_find( loc_hld{iS}(:,4), san_fra_loc );
    loc_ind{iS,4} = sum( [numel(loc_ind{iS,2}) numel(loc_ind{iS,3})] ) > 0;
    loc_ind{iS,5} = [ loc_ind{iS,2}' loc_ind{iS,3} ];
    
end

%%
out_tbl = cell( size(sbj_nme_use,1), 4 );
dub_tbl = cell( size(sbj_nme_use,1), 2 );

for iS = 1:size(sbj_nme_use,1)
    
    sbj_nme_ind = find( strcmpi( sbj_nme(:,1), sbj_nme_use{iS,1}) );
    
    out_tbl{iS,1} = sbj_nme_use{iS, 1};
    
    if numel(loc_ind{iS,5}) == 0
        out_tbl{iS,2}  = 'Empty';
        out_tbl{iS,3} = 'N/A';
        out_tbl{iS,4} = 'N/A';
    elseif numel(loc_ind{iS,5}) > 1
        
        dbl_ind = find( strcmpi( sbj_dbl_mri(:,1), sbj_nme_use{iS,1} ));
        
        if ~isempty(dbl_ind)
            chs_hld = string_find( loc_hld{iS}(:,1), sbj_dbl_mri{ dbl_ind, 2} );
            chs_loc = intersect( loc_ind{iS,5}, chs_hld);
            
            end_pnt = strfind( loc_hld{iS}{ chs_loc , 1 }, '/')-1;
            out_tbl{iS,2}  = loc_hld{iS}{ chs_loc , 1 }(1:end_pnt);
            out_tbl{iS,3} = loc_hld{iS}{ chs_loc , 2 };
            out_tbl{iS,4} = loc_hld{iS}{ chs_loc , 4 };
            
        else
            out_tbl{iS,14}  = 'Multiple';
            out_tbl{iS,15} = 'N/A';
            out_tbl{iS,16} = 'N/A';
            
            dub_tbl{iS,1} = sbj_nme_use{iS,1};
            dub_tbl{iS,2} = strcat( loc_hld{iS}( loc_ind{iS,5} , 1), '; ');
            dub_tbl{iS,2} = [dub_tbl{iS,2}{:}];
        end
        
        
    else
        end_pnt        = strfind( loc_hld{iS}{ loc_ind{iS,5} , 1 }, '/')-1;
        out_tbl{iS,2}  = loc_hld{iS}{ loc_ind{iS,5} , 1 }(1:end_pnt);
        out_tbl{iS,3} = loc_hld{iS}{ loc_ind{iS,5} , 2 };
        out_tbl{iS,4} = loc_hld{iS}{ loc_ind{iS,5} , 4 };
    end
    
end

dub_tbl( cellfun( @isempty, dub_tbl(:,1)), :) = [];
cell2csv([ prj_dir '/' 'Multiple_FSURFs_mri.csv'],dub_tbl);

mss_tbl = out_tbl;
mss_tbl( ~cellfun( @isempty, mss_tbl(:,2)), :) = [];

hav_tbl = out_tbl;
hav_tbl( cellfun( @isempty, hav_tbl(:,2)), :)   = [];
hav_tbl( strcmpi( hav_tbl(:,2),'Multiple'), :) = [];

prb_tbl = hav_tbl;
prb_tbl( ~strcmpi(prb_tbl(:,3), 'N/A'), :) = [];

hav_tbl( strcmpi(hav_tbl(:,3), 'N/A'), :) = [];
cell2csv([ prj_dir '/' 'FSURFs_mri.csv'],hav_tbl);

[ size(sbj_nme_use,1) size(hav_tbl,1) size(mss_tbl,1) size(prb_tbl,1) size(dub_tbl,1) ]

%% Create commands
frs_dir = '/space/syn09/1/data/MMILDB/MCD_RSI/raw/';
sec_dir = out_dir;

prn_cmd_mri = '';
for iS = 1:size(hav_tbl,1)
    
    prn_cmd_mri = sprintf('\n\n%s########################\n\nmkdir %s/%s/\n', prn_cmd_mri, sec_dir, hav_tbl{iS,1});
    prn_cmd_mri = sprintf('%smkdir %s/%s/mri/\n', prn_cmd_mri, sec_dir, hav_tbl{iS,1});
    prn_cmd_mri = sprintf('%scp -rL %s/%s/%s/*.dcm %s/%s/mri/\n', prn_cmd_mri, frs_dir, hav_tbl{iS,2}, hav_tbl{iS,3}, sec_dir, hav_tbl{iS,1});
    prn_cmd_mri = sprintf('%sls -l %s/%s/mri/ | wc -l\n\n\n\n', prn_cmd_mri, sec_dir, hav_tbl{iS,1});
    
end
cell2csv( [ out_dir '/' 'mri_copy_commands.csv'], {prn_cmd_mri});

%% Check Match between MRI & DTI
dti_hld = mmil_readtext([ prj_dir '/' 'FSURFs_dti.csv']);
mri_hld = mmil_readtext([ prj_dir '/' 'FSURFs_mri.csv']);

for iS = 1:size( dti_hld, 1)
    mri_loc = strcmpi( dti_hld{iS,1}, mri_hld(:,1));  
    
    if all(mri_loc==0)
        fprintf('\n\n#############\n%s: MRI MISSING!\n#############\n\n', dti_hld{iS,1})
    else
        mri_dti_mtc = strcmpi( dti_hld{iS,14}, mri_hld{ mri_loc, 2});
        
        if mri_dti_mtc
            fprintf('%s: Matches!\n', dti_hld{iS,1})
        else
            fprintf('\n\n#############\n%s: PROBLEM\n%s\n%s\n\n#############\n\n', dti_hld{iS,1}, dti_hld{iS,14}, mri_hld{ mri_loc, 2})
        end
    end
    
end






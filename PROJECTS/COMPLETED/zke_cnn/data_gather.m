clear; clc;

prj_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/CoAuthor/Zeke/MachineLearningAnnals';

%%
sbj_nme = mmil_readtext([prj_dir '/' 'Nonlesional_Engel_integrated.csv']);
sbj_add = mmil_readtext([prj_dir '/' 'MultimodalImagingOfL-3TDTIClinicalVariabl_DATA_LABELS_2020-05-20_1504.csv']);
sbj_dbl = mmil_readtext([prj_dir '/' 'Multiple_FSURFs_acm.csv']);

dta_loc = mmil_readtext([prj_dir '/' 'dti_search_noapostrophe.csv']);

for iS = 2:size(sbj_nme,1)
   
    row_ind = string_find( dta_loc(:,1),sbj_nme{iS} );
    
    loc_hld{iS} = cat( 1, dta_loc(row_ind, [1 2 4 7 ] ) );
    
end

[ sbj_nme(:,1) num2cell(cellfun( @(x) size(x,1), loc_hld ))' ];

% missing
% multiple

loc_ind = cell( size(sbj_nme,1), 5 );
for iS = 2:size(sbj_nme,1)
    
    san_die_loc = '30-dir DTI Pepolar #1';
    san_fra_loc = '(B=1000 DIR=30)';
    
    loc_ind{iS,1} = sbj_nme{iS,1};
    loc_ind{iS,2} = string_find( loc_hld{iS}(:,4), san_die_loc );
    loc_ind{iS,3} = string_find( loc_hld{iS}(:,4), san_fra_loc );
    loc_ind{iS,4} = sum( [numel(loc_ind{iS,2}) numel(loc_ind{iS,3})] ) > 0;
    loc_ind{iS,5} = [ loc_ind{iS,2} loc_ind{iS,3} ];
     
end

%%
out_tbl = cell( size(sbj_nme,1)-1, 11 );
dub_tbl = cell( size(sbj_nme,1)-1, 2 );
for iS = 2:size(sbj_nme,1)
    
    mtc_ind = find( strcmpi( sbj_add(:,1), sbj_nme{iS,1} ));
    
    out_tbl{iS-1,1} = sbj_nme{iS,1};
    
    if isempty(sbj_nme{iS,12})
        out_tbl{iS-1,2} = '';
    else
        out_tbl{iS-1,2} = sbj_nme{iS,12};
    end
    
    out_tbl{iS-1,3} = sbj_nme{iS,13};
        
    out_tbl{iS-1,4} = sbj_nme{iS,2};
    
    if sbj_nme{iS,3}==1
        out_tbl{iS-1,5} = 'MTS';
    else
        out_tbl{iS-1,5} = 'non-MTS';
    end
    
    if ~isempty(mtc_ind)
        out_tbl{iS-1,6} = sbj_add{ mtc_ind, 2 };
        out_tbl{iS-1,7} = sbj_add{ mtc_ind, 3 };
        out_tbl{iS-1,8} = sbj_add{ mtc_ind, 4 };
    end
    
    if numel(loc_ind{iS,5}) == 0
        out_tbl{iS-1,9} = 'Empty';    
        out_tbl{iS-1,10} = 'N/A';
        out_tbl{iS-1,11} = 'N/A';
    elseif numel(loc_ind{iS,5}) > 1
        
        dbl_ind = find( strcmpi( sbj_dbl(:,1), sbj_nme{iS,1} ));
        chs_hld = find( strcmpi( loc_hld{iS}(:,1), sbj_dbl{ dbl_ind, 2} ));
        chs_loc = intersect( loc_ind{iS,5}, chs_hld);
        
        end_pnt = strfind( loc_hld{iS}{ chs_loc , 1 }, '/')-1;
        out_tbl{iS-1,9}  = loc_hld{iS}{ chs_loc , 1 }(1:end_pnt);   
        out_tbl{iS-1,10} = loc_hld{iS}{ chs_loc , 2 };
        out_tbl{iS-1,11} = loc_hld{iS}{ chs_loc , 4 }; 
        
%         out_tbl{iS-1,9} = 'Multiple';   
%         out_tbl{iS-1,10} = 'N/A';
%         out_tbl{iS-1,11} = 'N/A';
%         
%         dub_tbl{iS-1,1} = sbj_nme{iS,1};
%         dub_tbl{iS-1,2} = strcat( loc_hld{iS}( loc_ind{iS,5} , 1), '; ');
%         dub_tbl{iS-1,2} = [dub_tbl{iS-1,2}{:}];
        
    else
        end_pnt = strfind( loc_hld{iS}{ loc_ind{iS,5} , 1 }, '/')-1;
        out_tbl{iS-1,9} = loc_hld{iS}{ loc_ind{iS,5} , 1 }(1:end_pnt);   
        out_tbl{iS-1,10} = loc_hld{iS}{ loc_ind{iS,5} , 2 };
        out_tbl{iS-1,11} = loc_hld{iS}{ loc_ind{iS,5} , 4 };   
    end
        
end

dub_tbl( cellfun( @isempty, dub_tbl(:,1)), :) = [];
cell2csv('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/CoAuthor/Zeke/Multiple_FSURFs.csv',dub_tbl);

mss_tbl = out_tbl;
mss_tbl( ~cellfun( @isempty, mss_tbl(:,2)), :) = [];

hav_tbl = out_tbl;
hav_tbl( cellfun( @isempty, hav_tbl(:,2)), :) = [];

prb_tbl = hav_tbl;
prb_tbl( ~strcmpi(prb_tbl(:,10), 'N/A'), :) = [];

hav_tbl( strcmpi(hav_tbl(:,10), 'N/A'), :) = [];

ttt = tabulate( [ hav_tbl(:,2)]); %; prb_tbl(:,2)]);
ttt = ttt(1:end,:);
[~, ttt_ind] = sort(ttt(:,1));
ttt = ttt(ttt_ind,:);

sum(cell2mat(ttt(1:4,2)))
sum(cell2mat(ttt(5:end,2)))

%% Creat csv
col_nme = { 'SubjID' 'Engel' 'Surgery' 'Seizure_Laterality' 'MTS_Status' 'Sex' 'Age_at_surgery' 'Age_at_seizure_onset' };
cell2csv( [ prj_dir '/' 'Subject_data.csv'], [ col_nme ; hav_tbl(:, 1:8) ]);

%% Create commands
frs_dir = '/space/syn09/1/data/MMILDB/MCD_RSI/raw/';
sec_dir = '/home/ekaestne/PROJECTS/DATA/ZEKE_SURGICAL_PROJECT';

for iS = 1:size(hav_tbl,1)

        prn_cmd = sprintf('\n\n########################\n\nmkdir %s/%s/\n', sec_dir, hav_tbl{iS,1});
        prn_cmd = sprintf('%scp -rL %s/%s/%s/*.dcm %s/%s/\n', prn_cmd, frs_dir, hav_tbl{iS,9}, hav_tbl{iS,10}, sec_dir, hav_tbl{iS,1});
        prn_cmd = sprintf('%sls -l %s/%s | wc -l\n\n\n\n', prn_cmd, sec_dir, hav_tbl{iS,1});
    
        fprintf('%s',prn_cmd)
        
end

%% Create commands, bvec/bval
frs_dir = '/space/syn09/1/data/MMILDB/MCD_RSI/proc_dti/';
sec_dir = '/home/ekaestne/PROJECTS/DATA/FSL_files';

for iS = 1:size(hav_tbl,1)

        prn_cmd = sprintf('\n\n########################\n\nmkdir %s/%s/\n', sec_dir, hav_tbl{iS,1});
        prn_cmd = sprintf('%scp -rL %s/%s/exportDTIforFSL/DTI1/*.txt %s/%s/\n', prn_cmd, frs_dir, strrep(hav_tbl{iS,9},'MRIRAW','DTIPROC'), sec_dir, hav_tbl{iS,1});
        prn_cmd = sprintf('%sls -l %s/%s | wc -l\n\n\n\n', prn_cmd, sec_dir, hav_tbl{iS,1});
    
        fprintf('%s',prn_cmd)
        
end


clear; clc;

dta_dir = '/space/mcdonald-syn01/1/data/MCD_MRI/EmorySLAH';

ejk_chk_dir('/home/ekaestne/Downloads/Emory_Slah_Explore/');

%%
sbj_dir = dir(dta_dir);
    sbj_dir = {sbj_dir(:).name};
    sbj_dir = sbj_dir(3:81);
  
%%    
sbj_tot_hld = cell(0);
scn_tot_hld = [];
for iS = 1:numel(sbj_dir)

    neu_dir = dir( [ dta_dir '/' sbj_dir{iS} ] );
    
    neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})) = {neu_dir(3:end).name};
    
    scn_tot_hld = [ scn_tot_hld numel(neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'}))) ];
    sbj_tot_hld = { sbj_tot_hld{:} neu_dir(3:end).name };
    
end

sbj_tot_hld = sort(sbj_tot_hld);
tbl_hld = tabulate(sbj_tot_hld);
sbj_tot_hld = sort(unique(sbj_tot_hld));

cell2csv( '/home/ekaestne/Downloads/Emory_Slah_Explore/emory_slah_table.csv', tbl_hld(:,1:2))

%%
tot_scn_num = max(scn_tot_hld);
tot_sbj_num = max(cell2mat(tbl_hld(:,2)));

scn_out_hld = cell( size(tbl_hld,1), tot_sbj_num+1);
    scn_out_hld(:,1) = tbl_hld(:,1);
sbj_out_hld = cell( numel(sbj_dir), tot_scn_num+1);
    sbj_out_hld(:,1) = sbj_dir';

for iS = 1:numel(sbj_dir)     
    for iSC = 1:numel(neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})))
        
    scn_num = dir( [ dta_dir '/' sbj_dir{iS} '/' neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){iSC} ] );
    scn_num = num2str(numel( string_find({scn_num(:).name},{'\.dcm$'}) ));
    
    % Subjects
    sbj_out_hld{iS,iSC+1} =[ neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){iSC} '; ' scn_num] ;
    
    % Scans
    scn_row = find(strcmpi( scn_out_hld(:,1), neu_hld.(mmil_spec_char(sbj_dir{iS},{'-' '.'})){iSC}));
    emp_col = find(cellfun(@isempty,scn_out_hld(scn_row,:)),1);
    scn_out_hld{ scn_row, emp_col} = [ sbj_dir{iS} ', ' scn_num];

    
    
    end
end

cell2csv( '/home/ekaestne/Downloads/Emory_Slah_Explore/emory_slah_scans.csv',    scn_out_hld)
cell2csv( '/home/ekaestne/Downloads/Emory_Slah_Explore/emory_slah_subjects.csv', sbj_out_hld)

%%




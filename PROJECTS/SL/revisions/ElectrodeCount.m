clear; clc;

prj_dat_hld = '/home/ekaestne/PROJECTS/OUTPUT/SL/';% '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/';
clr_fld     = [prj_dat_hld ];

sbj_nme_hld = mmil_readtext([clr_fld 'subjects']);
sbj_clr_hld = mmil_readtext([clr_fld 'subjects_clerical']);
sbj_dat_hld = mmil_readtext([clr_fld 'subjects_data']);

%% Run Subjects
out_dta = cell(numel(sbj_nme_hld),3);
for sbj_num = 1:numel(sbj_nme_hld);
   
   out_dta{sbj_num,1} = sbj_nme_hld{sbj_num};
    
   try
       lhs_ele = mmil_readtext([ prj_dat_hld '/' 'electrode_location_files' '/' sbj_nme_hld{sbj_num} '/' 'output' '/' sbj_nme_hld{sbj_num}  '_lhs_ecog' ]);
   catch
       lhs_ele = cell(0);
   end
   
   try
       rhs_ele = mmil_readtext([ prj_dat_hld '/' 'electrode_location_files' '/' sbj_nme_hld{sbj_num} '/' 'output' '/' sbj_nme_hld{sbj_num}  '_rhs_ecog' ]);
   catch
       rhs_ele = cell(0);
   end
   
   out_dta{sbj_num,2} = size(lhs_ele,1);
   out_dta{sbj_num,3} = size(rhs_ele,1);
   
end

sum(cell2mat(out_dta(:,2)))
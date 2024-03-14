sbj_nme_lng = mmil_readtext('/home/ekaestner/Dropbox/McDonald Lab/Erik/01_EndlessApps/2024/03_ML_R01/work/subjects.csv');
sbj_nme_lng = sbj_nme_lng(1:2:end,:);

fcfg = [];
fcfg.dta_loc = '/home/ekaestner/Dropbox/McDonald Lab/Erik/01_EndlessApps/2024/03_ML_R01/work/Diffusion_CNN_ClinicalData.csv';
[cln_dta,cln_sbj,cln_col] = ejk_dta_frm(fcfg);

%% Shorten
sbj_nme = cellfun(@(x) x(7),regexp(sbj_nme_lng,'/','split'));

%% Find outcomes
sbj_col = strcmpi(cln_col,'New ID');
cln_var = ismember(cln_col,{ 'Outcome' 'Binary Outcome' });

out_dta = [ sbj_nme_lng cell(numel(sbj_nme),sum(cln_var)) ];
for iS = 1:numel(sbj_nme)    
    sbj_row = strcmpi(cln_dta(:,sbj_col),sbj_nme{iS});
    if any(sbj_row)
        out_dta(iS,2:end) = cln_dta(sbj_row,cln_var);
    else
        fprintf('MISSING: %s \n',sbj_nme{iS})
    end
end

cell2csv('/home/ekaestner/Dropbox/McDonald Lab/Erik/01_EndlessApps/2024/03_ML_R01/work/Outcomes4Leo.csv',out_dta);

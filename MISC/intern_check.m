clear; clc;

dta_dir = '/home/ekaestner/Dropbox/McDonald Lab/RA_Stuff/HighTechHigh_Interns_2022';

int_nme = { 'Bianca'               'Solei' };
int_fle = { 'Bianca_Database.xlsx' 'Database_Solei.xlsx' };

sbj_num = 1:5;

bia_dta = readtable([ dta_dir '/' int_fle{1} ]);
    bia_dta_var = table2cell(bia_dta(:,1));
    bia_dta_sbj = bia_dta.Properties.VariableNames(2:end);
    bia_dta     = table2cell(bia_dta(:,2:end));
    
sol_dta = readtable([ dta_dir '/' int_fle{1} ]);
    sol_dta_var = table2cell(sol_dta(:,1));
    sol_dta_sbj = sol_dta.Properties.VariableNames(2:end);
    sol_dta     = table2cell(sol_dta(:,2:end));
    
%% Match/MisMatch
out_mtc = { [int_nme{1} '_' 'Subject']  [int_nme{2} '_' 'Subject']  ...
            [int_nme{1} '_' 'Variable'] [int_nme{2} '_' 'Variable'] ...
            [int_nme{1} '_' 'Value']    [int_nme{2} '_' 'Value'] };
for iS = 1:numel(sbj_num)
    
    for iV = 1:numel(bia_dta_var)
        
        bia_num = isnumeric(bia_dta{iV,sbj_num(iS)});
        sol_num = isnumeric(bia_dta{iV,sbj_num(iS)});
        
        if bia_num && sol_num
            bia_nan = isnan(bia_dta{iV,sbj_num(iS)});
            sol_nan = isnan(sol_dta{iV,sbj_num(iS)}); 
            
            if ~(bia_nan && sol_nan) && (bia_dta{iV,sbj_num(iS)} ~= sol_dta{iV,sbj_num(iS)})
                out_mtc = [ out_mtc ; bia_dta_sbj(sbj_num(iS)) sol_dta_sbj(sbj_num(iS)) bia_dta_var(iV) sol_dta_var(iV) bia_dta(iV,sbj_num(iS)) sol_dta(iV,sbj_num(iS)) ];
            end
        elseif  ~bia_num && ~sol_num
            if ~strcmpi(bia_dta{iV,sbj_num(iS)},sol_dta{iV,sbj_num(iS)})
               out_mtc = [ out_mtc ; bia_dta_sbj(sbj_num(iS)) sol_dta_sbj(sbj_num(iS)) bia_dta_var(iV) sol_dta_var(iV) bia_dta(iV,sbj_num(iS)) sol_dta(iV,sbj_num(iS)) ]; 
            end
        elseif  bia_num ~= sol_num
            out_mtc = [ out_mtc ; bia_dta_sbj(sbj_num(iS)) sol_dta_sbj(sbj_num(iS)) bia_dta_var(iV) sol_dta_var(iV) bia_dta(iV,sbj_num(iS)) sol_dta(iV,sbj_num(iS)) ];
        end        
    end    
end

cell2csv([ dta_dir '/' 'database_check' '/' 'Mismatch.csv' ],out_mtc)

%% MissingValues
for iS = 1:numel(sbj_num)
    out_mss = { [int_nme{1} '_' 'Subject']  [int_nme{2} '_' 'Subject']  ...
        [int_nme{1} '_' 'Variable'] [int_nme{2} '_' 'Variable'] ...
        [int_nme{1} '_' 'Value']    [int_nme{2} '_' 'Value'] };
    for iV = 1:numel(bia_dta_var)
        bia_num = isnumeric(bia_dta{iV,sbj_num(iS)});
        sol_num = isnumeric(bia_dta{iV,sbj_num(iS)});
        if bia_num && sol_num
            bia_nan = isnan(bia_dta{iV,sbj_num(iS)});
            sol_nan = isnan(sol_dta{iV,sbj_num(iS)});
            if bia_nan && sol_nan
                out_mss = [ out_mss ; bia_dta_sbj(sbj_num(iS)) sol_dta_sbj(sbj_num(iS)) bia_dta_var(iV) sol_dta_var(iV) bia_dta(iV,sbj_num(iS)) sol_dta(iV,sbj_num(iS)) ];
            elseif bia_nan ~= sol_nan
                error()
            end
        elseif  bia_num ~= sol_num
            error()
        end
    end
    cell2csv([ dta_dir '/' 'database_check' '/' 'MissingValues' '/' bia_dta_sbj{iS} '.csv' ],out_mss)
end

%% UnscoredValues
for iS = 1:numel(sbj_num)
    out_mss = { [int_nme{1} '_' 'Subject']  [int_nme{2} '_' 'Subject']  ...
        [int_nme{1} '_' 'Variable'] [int_nme{2} '_' 'Variable'] ...
        [int_nme{1} '_' 'Value']    [int_nme{2} '_' 'Value'] };
    for iV = 1:numel(bia_dta_var)
        bia_num = isnumeric(bia_dta{iV,sbj_num(iS)});
        sol_num = isnumeric(bia_dta{iV,sbj_num(iS)});
        if bia_num && sol_num
            bia_val = bia_dta{iV,sbj_num(iS)};
            sol_val = sol_dta{iV,sbj_num(iS)};
            if bia_val==-9999 && sol_val==-9999
                out_mss = [ out_mss ; bia_dta_sbj(sbj_num(iS)) sol_dta_sbj(sbj_num(iS)) bia_dta_var(iV) sol_dta_var(iV) bia_dta(iV,sbj_num(iS)) sol_dta(iV,sbj_num(iS)) ];
            end
        elseif  bia_num ~= sol_num
            error()
        end
    end
    cell2csv([ dta_dir '/' 'database_check' '/' 'UnscoredValues' '/' bia_dta_sbj{iS} '.csv' ],out_mss)
end






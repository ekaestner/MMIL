%% Grainsize 
load(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/figure4/grainsize' '/' 'tme_pnt.mat'])
tme_hld = 0:0.020:2;
for iH = 1:2
    for iR = 1:size(tme_plt{iH},1)
        for iC = 2:size(tme_plt{iH},2)
            if ~isempty(tme_plt{iH}{iR,iC})
                for iTM = 1:numel(tme_plt{iH}{iR,iC})
                    tme_plt{iH}{iR,iC}(iTM) = tme_hld(dsearchn(tme_hld',tme_plt{iH}{iR,iC}(iTM)));
                end
            end
        end
    end
end

fcfg.cmb_reg = { { 'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform'} ...
                 { 'lateraloccipital' } ...
                 { 'caudal-ITG' 'middle-ITG' 'rostral-ITG' } ...
                 { 'caudal-MTG' 'middle-MTG' 'rostral-MTG' } ...
                 { 'caudal-STG' 'middle-STG' 'rostral-STG' } ...
                 { 'inferior-precentral' 'middle-precentral' } ...
                 { 'parstriangularis' } ...
                 { 'parsopercularis' } };
fcfg.reg_nme = { 'Fusiform' ...
                 'Lateral Occipital' ...
                 'ITG' ...
                 'MTG' ...
                 'STG' ...
                 'Precentral' ...
                 'Pars Triangularis' ...
                 'Pars Operculatris' };
   
eff_lbl = { 'Letter' 'Word' 'False-Font' };
             
% stats %%%%%%%%%%%%%%%%
% put together regions
iH = 1;
for iC = 1:numel(fcfg.reg_nme)
    tab_hld{iC,1} = fcfg.reg_nme{iC};
    for iE = 2:size(tme_plt{iH},2)
        tab_hld{iC,iE} = [tme_plt{iH}{ismember(tme_plt{iH}(:,1),fcfg.cmb_reg{iC}),iE}];
        tab_hld{iC,iE} = tab_hld{iC,iE}(tab_hld{iC,iE}>0.100);
    end
end
tme_plt{iH} = tab_hld;

iH = 2;
for iC = 1:numel(fcfg.reg_nme)
    tab_hld{iC,1} = fcfg.reg_nme{iC};
    for iE = 2:size(tme_plt{iH},2)
        tab_hld{iC,iE} = [tme_plt{iH}{ismember(tme_plt{iH}(:,1),fcfg.cmb_reg{iC}),iE}];
    end
end
tme_plt{iH} = tab_hld;

% Interhemisphere
hms_tme = cell(size(tme_plt{1},1)+1,size(tme_plt,2));
for iR = 1:size(tme_plt{1},1)
    hms_tme{iR+1,1} = tme_plt{1}{iR,1};
    for iC = 2:size(tme_plt{1},2)
        hms_tme{1,iC} = eff_lbl{iC-1};
        if ~isempty(tme_plt{1}{iR,iC}) && ~isempty(tme_plt{2}{iR,iC})
        
            act_pvl = ranksum(tme_plt{1}{iR,iC},tme_plt{2}{iR,iC});
            bst_pvl = ranksum(ones(numel(tme_plt{1}{iR,iC}),1),ones(numel(tme_plt{2}{iR,iC}),1)*2);

            hms_tme{iR+1,iC} = [num2str(sprintf('%.2f',roundsd(act_pvl,2))) ' (' num2str(sprintf('%.2f',roundsd(bst_pvl,2))) ')'];
            
        else
            hms_tme{iR+1,iC} = '-';
        end        
    end
end

% Between Effects
cmb_col = {[1 2] [1 3] [2 3]};

for iH = 1:2
    cmb_tme{iH} = cell(size(tme_plt{iH},1)+1,size(cmb_col,2)+1);
    for iR = 1:size(tme_plt{1},1)
        cmb_tme{iH}{iR+1,1} = tme_plt{1}{iR,1};
        for iC = 1:size(cmb_col,2)
            cmb_tme{iH}{1,iC+1} = [eff_lbl{cmb_col{iC}(1)} '/' eff_lbl{cmb_col{iC}(2)} ];
            if ~isempty(tme_plt{iH}{iR,cmb_col{iC}(1)+1}) && ~isempty(tme_plt{iH}{iR,cmb_col{iC}(2)+1})
                
                act_pvl = ranksum(tme_plt{iH}{iR,cmb_col{iC}(1)+1},tme_plt{iH}{iR,cmb_col{iC}(2)+1});
                bst_pvl = ranksum(ones(numel(tme_plt{iH}{iR,cmb_col{iC}(1)+1}),1),ones(numel(tme_plt{iH}{iR,cmb_col{iC}(2)+1}),1)*2);
                
                cmb_tme{iH}{iR+1,iC+1} = [num2str(sprintf('%.2f',roundsd(act_pvl,2))) ' (' num2str(sprintf('%.2f',roundsd(bst_pvl,2))) ')'];
                
            else
                cmb_tme{iH}{iR+1,iC+1} = '-';
            end
        end
    end
end

% Between Regions
lhs_ltr = cell(size(tme_plt{iH},1),size(tme_plt{iH},1)+1);
lhs_wrd = cell(size(tme_plt{iH},1),size(tme_plt{iH},1)+1);
lhs_fls = cell(size(tme_plt{iH},1),size(tme_plt{iH},1)+1);

rhs_ltr = cell(size(tme_plt{iH},1),size(tme_plt{iH},1)+1);
rhs_wrd = cell(size(tme_plt{iH},1),size(tme_plt{iH},1)+1);
rhs_fls = cell(size(tme_plt{iH},1),size(tme_plt{iH},1)+1);

for iR = 2:size(tme_plt{iH},1)+1
    for iC = 2:size(tme_plt{iH},1)+1
        
        % lhs
        lhs_ltr{iR,1} = tme_plt{1}{iR-1,1};
        lhs_ltr{1,iC} = tme_plt{1}{iC-1,1};
        if ~isempty(tme_plt{1}{iR-1,2}) && ~isempty(tme_plt{1}{iC-1,2}) && iC ~= iR
            
            act_pvl = ranksum(tme_plt{1}{iR-1,2},tme_plt{1}{iC-1,2});
            bst_pvl = ranksum(ones(numel(tme_plt{1}{iR-1,2}),1),ones(numel(tme_plt{1}{iC-1,2}),1)*2);
            
            lhs_ltr{iR,iC} = [num2str(sprintf('%.2f',roundsd(act_pvl,2))) ' (' num2str(sprintf('%.2f',roundsd(bst_pvl,2))) ')'];
        else
            lhs_ltr{iR,iC} = '-';
        end
        
        lhs_wrd{iR,1} = tme_plt{1}{iR-1,1};
        lhs_wrd{1,iC} = tme_plt{1}{iC-1,1};
        if ~isempty(tme_plt{1}{iR-1,3}) && ~isempty(tme_plt{1}{iC-1,3}) && iC ~= iR
            
            act_pvl = ranksum(tme_plt{1}{iR-1,3},tme_plt{1}{iC-1,3});
            bst_pvl = ranksum(ones(numel(tme_plt{1}{iR-1,3}),1),ones(numel(tme_plt{1}{iC-1,3}),1)*2);
            
            lhs_wrd{iR,iC} = [num2str(sprintf('%.2f',roundsd(act_pvl,2))) ' (' num2str(sprintf('%.2f',roundsd(bst_pvl,2))) ')'];
        else
            lhs_wrd{iR,iC} = '-';
        end
        
        lhs_fls{iR,1} = tme_plt{1}{iR-1,1};
        lhs_fls{1,iC} = tme_plt{1}{iC-1,1};
        if ~isempty(tme_plt{1}{iR-1,4}) && ~isempty(tme_plt{1}{iC-1,4}) && iC ~= iR
            
            act_pvl = ranksum(tme_plt{1}{iR-1,4},tme_plt{1}{iC-1,4});
            bst_pvl = ranksum(ones(numel(tme_plt{1}{iR-1,4}),1),ones(numel(tme_plt{1}{iC-1,4}),1)*2);
            
            lhs_fls{iR,iC} = [num2str(sprintf('%.2f',roundsd(act_pvl,2))) ' (' num2str(sprintf('%.2f',roundsd(bst_pvl,2))) ')'];
        else
            lhs_fls{iR,iC} = '-';
        end
        
        % rhs
        rhs_ltr{iR,1} = tme_plt{2}{iR-1,1};
        rhs_ltr{1,iC} = tme_plt{2}{iC-1,1};
        if ~isempty(tme_plt{2}{iR-1,2}) && ~isempty(tme_plt{2}{iC-1,2}) && iC ~= iR
            
            act_pvl = ranksum(tme_plt{2}{iR-1,2},tme_plt{2}{iC-1,2});
            bst_pvl = ranksum(ones(numel(tme_plt{2}{iR-1,2}),1),ones(numel(tme_plt{2}{iC-1,2}),1)*2);
            
            rhs_ltr{iR,iC} = [num2str(sprintf('%.2f',roundsd(act_pvl,2))) ' (' num2str(sprintf('%.2f',roundsd(bst_pvl,2))) ')'];
        else
            rhs_ltr{iR,iC} = '-';
        end
        
        rhs_wrd{iR,1} = tme_plt{2}{iR-1,1};
        rhs_wrd{1,iC} = tme_plt{2}{iC-1,1};
        if ~isempty(tme_plt{2}{iR-1,3}) && ~isempty(tme_plt{2}{iC-1,3}) && iC ~= iR
            
            act_pvl = ranksum(tme_plt{2}{iR-1,3},tme_plt{2}{iC-1,3});
            bst_pvl = ranksum(ones(numel(tme_plt{2}{iR-1,3}),1),ones(numel(tme_plt{2}{iC-1,3}),1)*2);
            
            rhs_wrd{iR,iC} = [num2str(sprintf('%.2f',roundsd(act_pvl,2))) ' (' num2str(sprintf('%.2f',roundsd(bst_pvl,2))) ')'];
        else
            rhs_wrd{iR,iC} = '-';
        end
        
        rhs_fls{iR,1} = tme_plt{2}{iR-1,1};
        rhs_fls{1,iC} = tme_plt{2}{iC-1,1};
        if ~isempty(tme_plt{2}{iR-1,4}) && ~isempty(tme_plt{2}{iC-1,4}) && iC ~= iR
            
            act_pvl = ranksum(tme_plt{2}{iR-1,4},tme_plt{2}{iC-1,4});
            bst_pvl = ranksum(ones(numel(tme_plt{2}{iR-1,4}),1),ones(numel(tme_plt{2}{iC-1,4}),1)*2);
            
            rhs_fls{iR,iC} = [num2str(sprintf('%.2f',roundsd(act_pvl,2))) ' (' num2str(sprintf('%.2f',roundsd(bst_pvl,2))) ')'];
        else
            rhs_fls{iR,iC} = '-';
        end
        
    end
end

%% LEXICO-SEMANTIC TIMING
lex_tme = mmil_readtext();











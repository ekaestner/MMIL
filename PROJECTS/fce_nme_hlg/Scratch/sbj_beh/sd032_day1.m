 
%% SD032_DY1
clear; clc;

out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Ideas/fce_nme_stm/PreliminaryData/HalgLab_FN';

dta_hld = mmil_readtext('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Ideas/fce_nme_stm/PreliminaryData/HalgLab_FN/SD032_day1.csv');
    dta_hld_col = dta_hld(1,:);
    dta_hld     = dta_hld(2:end,:);

col_int = [ find(strcmpi(dta_hld_col,'list'))  find(strcmpi(dta_hld_col,'accuracy')) find(strcmpi(dta_hld_col,'RT'))];
    
trl_num_typ = {'Rep0' 'Rep1' 'Rep2' 'Rep3'};

lst_nme = unique(dta_hld(:,strcmpi(dta_hld_col,'list')));

out_csv = [];
for iL = 1:numel(lst_nme)
    lst_ind = find(strcmpi(dta_hld(:,strcmpi(dta_hld_col,'list')),lst_nme{iL}));
    
    op0_ind = lst_ind(find(strcmpi(dta_hld(lst_ind,strcmpi(dta_hld_col,'condition')),'io1')));
        num_op0_ind = numel(op0_ind);
        num_op0_brk = (num_op0_ind/3)*2;
    op1_ind = op0_ind(num_op0_brk+1:num_op0_ind);
    op0_ind = [op0_ind(3:2:num_op0_brk) ; op0_ind(num_op0_brk)];
    
    op2_ind = lst_ind(strcmpi(dta_hld(lst_ind,strcmpi(dta_hld_col,'condition')),'o2'));
    op3_ind = lst_ind(strcmpi(dta_hld(lst_ind,strcmpi(dta_hld_col,'condition')),'o3'));
    
    out_csv = [out_csv ; ...
               repmat({trl_num_typ{1}},numel(op0_ind),1) dta_hld(op0_ind,col_int) ; ...
               repmat({trl_num_typ{2}},numel(op1_ind),1) dta_hld(op1_ind,col_int) ; ...
               repmat({trl_num_typ{3}},numel(op2_ind),1) dta_hld(op2_ind,col_int) ; ...
               repmat({trl_num_typ{4}},numel(op3_ind),1) dta_hld(op3_ind,col_int) ];
    
end

cell2csv( [out_dir '/' 'SD032_DY1.csv'], out_csv);
clear out_csv

%% SD030_DY1
clear; clc;

out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Ideas/fce_nme_stm/PreliminaryData/HalgLab_FN';

dta_hld = mmil_readtext('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Ideas/fce_nme_stm/PreliminaryData/HalgLab_FN/SD030.csv');
    dta_hld_col = dta_hld(1,:);
    dta_hld     = dta_hld(2:end,:);

col_int = [ find(strcmpi(dta_hld_col,'list'))  find(strcmpi(dta_hld_col,'accuracy')) find(strcmpi(dta_hld_col,'RT'))];
    
trl_num_typ = {'Rep0' 'Rep1' 'Rep2' 'Rep3'};

lst_nme = unique(dta_hld(:,strcmpi(dta_hld_col,'list')));

out_csv = [];
for iL = 1:numel(lst_nme)
    lst_ind = find(strcmpi(dta_hld(:,strcmpi(dta_hld_col,'list')),lst_nme{iL}));
    
    op0_ind = lst_ind(find(strcmpi(dta_hld(lst_ind,strcmpi(dta_hld_col,'condition')),'io1')));
        num_op0_ind = numel(op0_ind);
        num_op0_brk = (num_op0_ind/3)*2;
    op1_ind = op0_ind(num_op0_brk+1:num_op0_ind);
    op0_ind = [op0_ind(3:2:num_op0_brk) ; op0_ind(num_op0_brk)];
    
    op2_ind = lst_ind(strcmpi(dta_hld(lst_ind,strcmpi(dta_hld_col,'condition')),'o2'));
    op3_ind = lst_ind(strcmpi(dta_hld(lst_ind,strcmpi(dta_hld_col,'condition')),'o3'));
    
    out_csv = [out_csv ; ...
               repmat({trl_num_typ{1}},numel(op0_ind),1) dta_hld(op0_ind,col_int) ; ...
               repmat({trl_num_typ{2}},numel(op1_ind),1) dta_hld(op1_ind,col_int) ; ...
               repmat({trl_num_typ{3}},numel(op2_ind),1) dta_hld(op2_ind,col_int) ; ...
               repmat({trl_num_typ{4}},numel(op3_ind),1) dta_hld(op3_ind,col_int) ];
    
end

cell2csv( [out_dir '/' 'SD030_DY1.csv'], out_csv);
clear out_csv

%% SD029_DY1
clear; clc;

out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Ideas/fce_nme_stm/PreliminaryData/HalgLab_FN';

dta_hld = mmil_readtext('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Ideas/fce_nme_stm/PreliminaryData/HalgLab_FN/SD029.csv');
    dta_hld_col = dta_hld(1,:);
    dta_hld     = dta_hld(2:end,:);

col_int = [ find(strcmpi(dta_hld_col,'list'))  find(strcmpi(dta_hld_col,'accuracy')) find(strcmpi(dta_hld_col,'RT'))];
    
trl_num_typ = {'Rep0' 'Rep1' 'Rep2' 'Rep3'};

lst_nme = unique(dta_hld(:,strcmpi(dta_hld_col,'list')));

out_csv = [];
for iL = 1:numel(lst_nme)
    lst_ind = find(strcmpi(dta_hld(:,strcmpi(dta_hld_col,'list')),lst_nme{iL}));
    
    op0_ind = lst_ind(find(strcmpi(dta_hld(lst_ind,strcmpi(dta_hld_col,'condition')),'io1')));
        num_op0_ind = numel(op0_ind);
        num_op0_brk = (num_op0_ind/3)*2;
    op1_ind = op0_ind(num_op0_brk+1:num_op0_ind);
    op0_ind = [op0_ind(3:2:num_op0_brk) ; op0_ind(num_op0_brk)];
    
    op2_ind = lst_ind(strcmpi(dta_hld(lst_ind,strcmpi(dta_hld_col,'condition')),'o2'));
    op3_ind = lst_ind(strcmpi(dta_hld(lst_ind,strcmpi(dta_hld_col,'condition')),'o3'));
    
    out_csv = [out_csv ; ...
               repmat({trl_num_typ{1}},numel(op0_ind),1) dta_hld(op0_ind,col_int) ; ...
               repmat({trl_num_typ{2}},numel(op1_ind),1) dta_hld(op1_ind,col_int) ; ...
               repmat({trl_num_typ{3}},numel(op2_ind),1) dta_hld(op2_ind,col_int) ; ...
               repmat({trl_num_typ{4}},numel(op3_ind),1) dta_hld(op3_ind,col_int) ];
    
end

cell2csv( [out_dir '/' 'SD029_DY1.csv'], out_csv);
clear out_csv

%% SD028_DY1
clear; clc;

out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Ideas/fce_nme_stm/PreliminaryData/HalgLab_FN';

dta_hld = mmil_readtext('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Ideas/fce_nme_stm/PreliminaryData/HalgLab_FN/SD028_day1.csv');
    dta_hld_col = dta_hld(1,:);
    dta_hld     = dta_hld(2:end,:);

col_int = [ find(strcmpi(dta_hld_col,'list'))  find(strcmpi(dta_hld_col,'accuracy')) find(strcmpi(dta_hld_col,'RT'))];
    
trl_num_typ = {'Rep0' 'Rep1' 'Rep2' 'Rep3'};

lst_nme = unique(dta_hld(:,strcmpi(dta_hld_col,'list')));

out_csv = [];
for iL = 1:numel(lst_nme)
    lst_ind = find(strcmpi(dta_hld(:,strcmpi(dta_hld_col,'list')),lst_nme{iL}));
    
    op0_ind = lst_ind(find(strcmpi(dta_hld(lst_ind,strcmpi(dta_hld_col,'condition')),'io1')));
        num_op0_ind = numel(op0_ind);
        num_op0_brk = (num_op0_ind/3)*2;
    op1_ind = op0_ind(num_op0_brk+1:num_op0_ind);
    op0_ind = [op0_ind(3:2:num_op0_brk) ; op0_ind(num_op0_brk)];
    
    op2_ind = lst_ind(strcmpi(dta_hld(lst_ind,strcmpi(dta_hld_col,'condition')),'o2'));
    op3_ind = lst_ind(strcmpi(dta_hld(lst_ind,strcmpi(dta_hld_col,'condition')),'o3'));
    
    out_csv = [out_csv ; ...
               repmat({trl_num_typ{1}},numel(op0_ind),1) dta_hld(op0_ind,col_int) ; ...
               repmat({trl_num_typ{2}},numel(op1_ind),1) dta_hld(op1_ind,col_int) ; ...
               repmat({trl_num_typ{3}},numel(op2_ind),1) dta_hld(op2_ind,col_int) ; ...
               repmat({trl_num_typ{4}},numel(op3_ind),1) dta_hld(op3_ind,col_int) ];
    
end

cell2csv( [out_dir '/' 'SD028_DY1.csv'], out_csv);
clear out_csv

%% OH025_DY1
clear; clc;

out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Ideas/fce_nme_stm/PreliminaryData/HalgLab_FN';

dta_hld = mmil_readtext('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Ideas/fce_nme_stm/PreliminaryData/HalgLab_FN/OHSU_HN_025.csv');
    dta_hld(162:end,:) = [];
    dta_hld_col = dta_hld(1,:);
    dta_hld     = dta_hld(2:end,:);

col_int = [ find(strcmpi(dta_hld_col,'list'))  find(strcmpi(dta_hld_col,'accuracy')) find(strcmpi(dta_hld_col,'RT'))];
    
trl_num_typ = {'Rep0' 'Rep1' 'Rep2' 'Rep3'};

lst_nme = unique(dta_hld(:,strcmpi(dta_hld_col,'list')));

out_csv = [];
for iL = 1:numel(lst_nme)
    lst_ind = find(strcmpi(dta_hld(:,strcmpi(dta_hld_col,'list')),lst_nme{iL}));
    
    op0_ind = lst_ind(find(strcmpi(dta_hld(lst_ind,strcmpi(dta_hld_col,'condition')),'io1')));
        num_op0_ind = numel(op0_ind);
        num_op0_brk = (num_op0_ind/3)*2;
    op1_ind = op0_ind(num_op0_brk+1:num_op0_ind);
    op0_ind = [op0_ind(3:2:num_op0_brk) ; op0_ind(num_op0_brk)];
    
    op2_ind = lst_ind(strcmpi(dta_hld(lst_ind,strcmpi(dta_hld_col,'condition')),'o2'));
    op3_ind = lst_ind(strcmpi(dta_hld(lst_ind,strcmpi(dta_hld_col,'condition')),'o3'));
    
    out_csv = [out_csv ; ...
               repmat({trl_num_typ{1}},numel(op0_ind),1) dta_hld(op0_ind,col_int) ; ...
               repmat({trl_num_typ{2}},numel(op1_ind),1) dta_hld(op1_ind,col_int) ; ...
               repmat({trl_num_typ{3}},numel(op2_ind),1) dta_hld(op2_ind,col_int) ; ...
               repmat({trl_num_typ{4}},numel(op3_ind),1) dta_hld(op3_ind,col_int) ];
    
end

cell2csv( [out_dir '/' 'OH025_DY1.csv'], out_csv);
clear out_csv

%% OH024_DY1
clear; clc;

out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Ideas/fce_nme_stm/PreliminaryData/HalgLab_FN';

dta_hld = mmil_readtext('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Ideas/fce_nme_stm/PreliminaryData/HalgLab_FN/OHSU_HN_024.csv');
    dta_hld(42:end,:) = [];
    dta_hld_col = dta_hld(1,:);
    dta_hld     = dta_hld(2:end,:);

col_int = [ find(strcmpi(dta_hld_col,'list'))  find(strcmpi(dta_hld_col,'accuracy')) find(strcmpi(dta_hld_col,'RT'))];
    
trl_num_typ = {'Rep0' 'Rep1' 'Rep2' 'Rep3'};

lst_nme = unique(dta_hld(:,strcmpi(dta_hld_col,'list')));

out_csv = [];
for iL = 1:numel(lst_nme)
    lst_ind = find(strcmpi(dta_hld(:,strcmpi(dta_hld_col,'list')),lst_nme{iL}));
    
    op0_ind = lst_ind(find(strcmpi(dta_hld(lst_ind,strcmpi(dta_hld_col,'condition')),'io1')));
        num_op0_ind = numel(op0_ind);
        num_op0_brk = (num_op0_ind/3)*2;
    op1_ind = op0_ind(num_op0_brk+1:num_op0_ind);
    op0_ind = [op0_ind(3:2:num_op0_brk) ; op0_ind(num_op0_brk)];
    
    op2_ind = lst_ind(strcmpi(dta_hld(lst_ind,strcmpi(dta_hld_col,'condition')),'o2'));
    op3_ind = lst_ind(strcmpi(dta_hld(lst_ind,strcmpi(dta_hld_col,'condition')),'o3'));
    
    out_csv = [out_csv ; ...
               repmat({trl_num_typ{1}},numel(op0_ind),1) dta_hld(op0_ind,col_int) ; ...
               repmat({trl_num_typ{2}},numel(op1_ind),1) dta_hld(op1_ind,col_int) ; ...
               repmat({trl_num_typ{3}},numel(op2_ind),1) dta_hld(op2_ind,col_int) ; ...
               repmat({trl_num_typ{4}},numel(op3_ind),1) dta_hld(op3_ind,col_int) ];
    
end

cell2csv( [out_dir '/' 'OH024_DY1.csv'], out_csv);
clear out_csv

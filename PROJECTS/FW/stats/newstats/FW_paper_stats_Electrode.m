clear; clc

% SETUP REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reg.oct_reg = { 'Caudal ITG' 'Middle ITG' 'Rostral ITG' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Lateral Occipital' };
reg.par_reg = { 'Supramarginal' 'Inferior Parietal' 'Superior Parietal' };
reg.tmp_reg = { 'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' };
reg.rol_reg = { 'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' };
reg.frn_reg = { 'rostral-middlefrontal' 'middle-middlefrontal' 'Caudal Middle Frontal' 'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis'};
reg.ifg_reg = { 'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis'};

cmb_reg.occ = { 'Lateral Occipital' };
cmb_reg.fus = { 'Caudal Fusiform' 'Middle Fusiform' };
cmb_reg.itg = { 'Caudal ITG' 'Middle ITG' };
cmb_reg.mtg = { 'Caudal MTG' 'Middle MTG' };
cmb_reg.stg = { 'Caudal STG' 'Middle STG' };
cmb_reg.par = { 'Inferior Parietal' 'Superior Parietal' };
cmb_reg.sup = { 'Supramarginal' };
cmb_reg.pre = { 'Inferior Precentral' 'Middle Precentral' };
cmb_reg.pre = { 'Inferior Precentral' 'Middle Precentral' };
cmb_reg.opc = { 'Pars Opercularis' };
cmb_reg.tri = { 'Pars Triangularis' };
% cmb_reg.orb = { 'Pars Orbitalis' };
cmb_reg.mid = { 'rostral-middlefrontal' 'middle-middlefrontal' 'Caudal Middle Frontal' };

%% Orth Overlap
nme = fieldnames(reg);
cmb_nme = fieldnames(cmb_reg);

tst_nme = 'pap_wrd_600';
hms     = { 'lhs' 'rhs' };
col     = { [2 3] [12 13] [22 23] };

tst.lhs = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/' tst_nme '/total/total_' tst_nme '_lhs_table_plot']);
tst.lhs = tst.lhs(2:end-1,:);
tmp_hld = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/' 'pap_con_600' '/total/total_' 'pap_con_600' '_lhs_table_plot']);
tst.lhs = [tst.lhs tmp_hld(2:end-1,[2 3])];

tst.rhs = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/' tst_nme '/total/total_' tst_nme '_rhs_table_plot']);
tst.rhs = tst.rhs(2:end-1,:);
tmp_hld = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/' 'pap_con_600' '/total/total_' 'pap_con_600' '_rhs_table_plot']);
tst.rhs = [tst.rhs tmp_hld(2:end-1,[2 3])];

% Pooled
reg_hld = cell(1,1);
for iG = 1:numel(col)
    for iR = 1:numel(nme)
        
        reg_1st = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], reg.(nme{iR}),'uni',0)));
        reg_2nd = find(ismember(tst.(hms{2})(:,1),cellfun(@(x) [hms{2} '_' x], reg.(nme{iR}),'uni',0)));
        
        ele_num_1st = 0; tot_num_1st = 0; ele_num_2nd = 0; tot_num_2nd = 0;
        for iN = 1:numel(reg_1st)
            ele_num_1st = ele_num_1st + tst.(hms{1}){reg_1st(iN),col{iG}(1)};
            tot_num_1st = tot_num_1st + tst.(hms{1}){reg_1st(iN),col{iG}(2)};
            
            ele_num_2nd = ele_num_2nd + tst.(hms{2}){reg_2nd(iN),col{iG}(1)};
            tot_num_2nd = tot_num_2nd + tst.(hms{2}){reg_2nd(iN),col{iG}(2)};
        end
        
        yyy = [ ones(1,ele_num_1st)   ones(1,tot_num_1st-ele_num_1st)   ...
            ones(1,ele_num_2nd)*2 ones(1,tot_num_2nd-ele_num_2nd)*2 ];
        
        zzz = [ ones(1,ele_num_1st) ones(1,tot_num_1st-ele_num_1st)*2 ...
            ones(1,ele_num_2nd) ones(1,tot_num_2nd-ele_num_2nd)*2 ];
        
        [sig,pvl] = FisherExactTest(yyy,zzz);
        
        reg_hld{iG}{iR,1} = nme{iR};
        reg_hld{iG}{iR,2} = sig;
        reg_hld{iG}{iR,3} = pvl;
        reg_hld{iG}{iR,4} = crosstab(yyy,zzz);
        if size(reg_hld{iG}{iR,4},2)==1; reg_hld{iG}{iR,4} = [0 reg_hld{iG}{iR,4}(1) ; 0 reg_hld{iG}{iR,4}(2)]; end
        
    end
end

pct_hld = cat(3,reg_hld{1}{:,4});
pct_bar = [squeeze(pct_hld(1,1,:))./(squeeze(pct_hld(1,2,:))+squeeze(pct_hld(1,1,:))) squeeze(pct_hld(2,1,:))./(squeeze(pct_hld(2,2,:))+squeeze(pct_hld(2,1,:)))];
figure()
bar(pct_bar)
set(gca,'XTickLabel',nme)

pct_hld = cat(3,reg_hld{2}{:,4});
pct_bar = [squeeze(pct_hld(1,1,:))./(squeeze(pct_hld(1,2,:))+squeeze(pct_hld(1,1,:))) squeeze(pct_hld(2,1,:))./(squeeze(pct_hld(2,2,:))+squeeze(pct_hld(2,1,:)))];
figure()
bar(pct_bar)
set(gca,'XTickLabel',nme)

pct_hld = cat(3,reg_hld{3}{:,4});
pct_bar = [squeeze(pct_hld(1,1,:))./(squeeze(pct_hld(1,2,:))+squeeze(pct_hld(1,1,:))) squeeze(pct_hld(2,1,:))./(squeeze(pct_hld(2,2,:))+squeeze(pct_hld(2,1,:)))];
figure()
bar(pct_bar)
set(gca,'XTickLabel',nme)

% Individual regions
cmb_hld = cell(1,1);
for iG = 1:numel(col)
    for iR = 1:numel(cmb_nme)
        
        reg_1st = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], cmb_reg.(cmb_nme{iR}),'uni',0)));
        reg_2nd = find(ismember(tst.(hms{2})(:,1),cellfun(@(x) [hms{2} '_' x], cmb_reg.(cmb_nme{iR}),'uni',0)));
        
        ele_num_1st = 0; tot_num_1st = 0; ele_num_2nd = 0; tot_num_2nd = 0;
        for iN = 1:numel(reg_1st)
            ele_num_1st = ele_num_1st + tst.(hms{1}){reg_1st(iN),col{iG}(1)};
            tot_num_1st = tot_num_1st + tst.(hms{1}){reg_1st(iN),col{iG}(2)};
            
            ele_num_2nd = ele_num_2nd + tst.(hms{2}){reg_2nd(iN),col{iG}(1)};
            tot_num_2nd = tot_num_2nd + tst.(hms{2}){reg_2nd(iN),col{iG}(2)};
        end
        
        yyy = [ ones(1,ele_num_1st)   ones(1,tot_num_1st-ele_num_1st)   ...
            ones(1,ele_num_2nd)*2 ones(1,tot_num_2nd-ele_num_2nd)*2 ];
        
        zzz = [ ones(1,ele_num_1st) ones(1,tot_num_1st-ele_num_1st)*2 ...
            ones(1,ele_num_2nd) ones(1,tot_num_2nd-ele_num_2nd)*2 ];
        
        [sig,pvl] = FisherExactTest(yyy,zzz);
        
        cmb_hld{iG}{iR,1} = cmb_nme{iR};
        cmb_hld{iG}{iR,2} = sig;
        cmb_hld{iG}{iR,3} = pvl;
        cmb_hld{iG}{iR,4} = crosstab(yyy,zzz);
        if size(cmb_hld{iG}{iR,4},2)==1; cmb_hld{iG}{iR,4} = [0 cmb_hld{iG}{iR,4}(1) ; 0 cmb_hld{iG}{iR,4}(2)]; end
        
    end
end

pct_hld = cat(3,cmb_hld{1}{:,4});
pct_bar = [squeeze(pct_hld(1,1,:))./(squeeze(pct_hld(1,2,:))+squeeze(pct_hld(1,1,:))) squeeze(pct_hld(2,1,:))./(squeeze(pct_hld(2,2,:))+squeeze(pct_hld(2,1,:)))];
figure()
bar(pct_bar)
title('Letter')
set(gca,'XTickLabel',cmb_nme)
'Letter';
cell2csv(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/tables' '/' 'ort_tbl' '/' 'ort_letter_num.csv'],[strcat('lhs_',cmb_hld{1}(:,1)) num2cell(squeeze(pct_hld(1,1,:))) strcat('rhs_',cmb_hld{1}(:,1)) num2cell(squeeze(pct_hld(2,1,:)))]);
cell2csv(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/tables' '/' 'ort_tbl' '/' 'ort_letter_sig.csv'],cmb_hld{1}(:,1:3));

pct_hld = cat(3,cmb_hld{2}{:,4});
pct_bar = [squeeze(pct_hld(1,1,:))./(squeeze(pct_hld(1,2,:))+squeeze(pct_hld(1,1,:))) squeeze(pct_hld(2,1,:))./(squeeze(pct_hld(2,2,:))+squeeze(pct_hld(2,1,:)))];
figure()
bar(pct_bar)
title('Word')
set(gca,'XTickLabel',cmb_nme)
'Word-Font';
cell2csv(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/tables' '/' 'ort_tbl' '/' 'ort_word_num.csv'],[strcat('lhs_',cmb_hld{2}(:,1)) num2cell(squeeze(pct_hld(1,1,:))) strcat('rhs_',cmb_hld{2}(:,1)) num2cell(squeeze(pct_hld(2,1,:)))])
cell2csv(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/tables' '/' 'ort_tbl' '/' 'ort_word_sig.csv'],cmb_hld{2}(:,1:3))

pct_hld = cat(3,cmb_hld{3}{:,4});
pct_bar = [squeeze(pct_hld(1,1,:))./(squeeze(pct_hld(1,2,:))+squeeze(pct_hld(1,1,:))) squeeze(pct_hld(2,1,:))./(squeeze(pct_hld(2,2,:))+squeeze(pct_hld(2,1,:)))];
figure()
bar(pct_bar)
title('False-Font')
set(gca,'XTickLabel',cmb_nme)
'False-Font';
cell2csv(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/tables' '/' 'ort_tbl' '/' 'ort_falsefont_num.csv'],[strcat('lhs_',cmb_hld{3}(:,1)) num2cell(squeeze(pct_hld(1,1,:))) strcat('rhs_',cmb_hld{3}(:,1)) num2cell(squeeze(pct_hld(2,1,:)))])
cell2csv(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/tables' '/' 'ort_tbl' '/' 'ort_falsefont_sig.csv'],cmb_hld{3}(:,1:3))

% COMPARING Within %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col     = { [2 3] [12 13] [22 23] };
clear rrr ccc crcrcrc

% Pooled
for iG = 1:numel(col)
    
    pll_lef_sig{iG} = cell(1,1);
    pll_lef_pvl{iG} = cell(1,1);
    pll_lef_crs{iG} = cell(1,1);
    ccc{iG}     = cell(0,0);
    rrr{iG}     = cell(0,0);
    crcrcrc{iG} = cell(0,0);
    
    for iR = 1:numel(nme)
        
        for iC = 1:numel(nme)-1
            
            if iR ~= iC
                reg_1st = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], reg.(nme{iR}),'uni',0)));
                reg_2nd = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], reg.(nme{iC}),'uni',0)));
                
                ele_num_1st = 0; tot_num_1st = 0; ele_num_2nd = 0; tot_num_2nd = 0;
                for iN = 1:numel(reg_1st)
                    ele_num_1st = ele_num_1st + tst.(hms{1}){reg_1st(iN),col{iG}(1)};
                    tot_num_1st = tot_num_1st + tst.(hms{1}){reg_1st(iN),col{iG}(2)};
                end
                for iN = 1:numel(reg_2nd)
                    ele_num_2nd = ele_num_2nd + tst.(hms{1}){reg_2nd(iN),col{iG}(1)};
                    tot_num_2nd = tot_num_2nd + tst.(hms{1}){reg_2nd(iN),col{iG}(2)};
                end
                
                yyy = [ ones(1,ele_num_1st)   ones(1,tot_num_1st-ele_num_1st)   ...
                    ones(1,ele_num_2nd)*2 ones(1,tot_num_2nd-ele_num_2nd)*2 ];
                
                zzz = [ ones(1,ele_num_1st) ones(1,tot_num_1st-ele_num_1st)*2 ...
                    ones(1,ele_num_2nd) ones(1,tot_num_2nd-ele_num_2nd)*2 ];
                
                [sig,pvl] = FisherExactTest(yyy,zzz);
                
                pll_lef_sig{iG}{iR+1,1} = nme{iR};
                pll_lef_sig{iG}{1,iC+1} = nme{iC};
                
                pll_lef_pvl{iG}{iR+1,1} = nme{iR};
                pll_lef_pvl{iG}{1,iC+1} = nme{iC};
                
                pll_lef_crs{iG}{iR+1,1} = nme{iR};
                pll_lef_crs{iG}{1,iC+1} = nme{iC};
                
                pll_lef_sig{iG}{iR+1,iC+1} = sig;
                pll_lef_pvl{iG}{iR+1,iC+1} = pvl;
                pll_lef_crs{iG}{iR+1,iC+1} = crosstab(yyy,zzz);
                if size(pll_lef_crs{iG}{iR+1,iC+1},2)==1; pll_lef_crs{iG}{iR+1,iC+1} = [0 pll_lef_crs{iG}{iR+1,iC+1}(1) ; 0 pll_lef_crs{iG}{iR+1,iC+1}(2)]; end      
                
                if ~sig && pvl>0.2 && ~any(strcmpi([nme{iR} nme{iC}],crcrcrc{iG})) &&  ~any(strcmpi([nme{iC} nme{iR}],crcrcrc{iG}))
                    rrr{iG}{end+1} = nme{iR};
                    ccc{iG}{end+1} = nme{iC};
                    crcrcrc{iG}{end+1} = [nme{iR} nme{iC}];
                end
                
            else
                pll_lef_crs{iG}{iR+1,iC+1} = '-';
            end
            
        end
        
        crs_hld         = crosstab(yyy,zzz);
        if size(crs_hld,2)==1; crs_hld = [0 crs_hld(1) ; 0 crs_hld(2)]; end
        pll_lef_bar{iG}(iR) = crs_hld(1,1)/crs_hld(1,2);
        
    end
    
    figure(); aaa{iG} = graph(rrr{iG},ccc{iG}); plot(aaa{iG});
    
end

figure()
bar(pll_lef_bar{1})
title('Letter')
set(gca,'XTickLabel',nme)

figure()
bar(pll_lef_bar{2})
title('Word')
set(gca,'XTickLabel',nme)

figure()
bar(pll_lef_bar{3})
title('False-Font')
set(gca,'XTickLabel',nme)

% Specific
clear rrr ccc crcrcrc

for iG = 1:numel(col)
    
    lef_sig{iG} = cell(1,1);
    lef_pvl{iG} = cell(1,1);
    lef_crs{iG} = cell(1,1);
    ccc{iG}     = cell(0,0);
    rrr{iG}     = cell(0,0);
    crcrcrc{iG} = cell(0,0);
    
    for iR = 1:numel(cmb_nme)
        
        for iC = 1:numel(cmb_nme)-1
            
            if iR ~= iC
                reg_1st = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], cmb_reg.(cmb_nme{iR}),'uni',0)));
                reg_2nd = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], cmb_reg.(cmb_nme{iC}),'uni',0)));
                
                ele_num_1st = 0; tot_num_1st = 0; ele_num_2nd = 0; tot_num_2nd = 0;
                for iN = 1:numel(reg_1st)
                    ele_num_1st = ele_num_1st + tst.(hms{1}){reg_1st(iN),col{iG}(1)};
                    tot_num_1st = tot_num_1st + tst.(hms{1}){reg_1st(iN),col{iG}(2)};
                end
                for iN = 1:numel(reg_2nd)
                    ele_num_2nd = ele_num_2nd + tst.(hms{1}){reg_2nd(iN),col{iG}(1)};
                    tot_num_2nd = tot_num_2nd + tst.(hms{1}){reg_2nd(iN),col{iG}(2)};
                end
                
                yyy = [ ones(1,ele_num_1st)   ones(1,tot_num_1st-ele_num_1st)   ...
                    ones(1,ele_num_2nd)*2 ones(1,tot_num_2nd-ele_num_2nd)*2 ];
                
                zzz = [ ones(1,ele_num_1st) ones(1,tot_num_1st-ele_num_1st)*2 ...
                    ones(1,ele_num_2nd) ones(1,tot_num_2nd-ele_num_2nd)*2 ];
                
                [sig,pvl] = FisherExactTest(yyy,zzz);
                
                lef_sig{iG}{iR+1,1} = cmb_nme{iR};
                lef_sig{iG}{1,iC+1} = cmb_nme{iC};
                
                lef_pvl{iG}{iR+1,1} = cmb_nme{iR};
                lef_pvl{iG}{1,iC+1} = cmb_nme{iC};
                
                lef_crs{iG}{iR+1,1} = cmb_nme{iR};
                lef_crs{iG}{1,iC+1} = cmb_nme{iC};
                
                lef_sig{iG}{iR+1,iC+1} = sig;
                lef_pvl{iG}{iR+1,iC+1} = pvl;
                lef_crs{iG}{iR+1,iC+1} = crosstab(yyy,zzz);
                if size(lef_crs{iG}{iR+1,iC+1},2)==1; lef_crs{iG}{iR+1,iC+1} = [0 lef_crs{iG}{iR+1,iC+1}(1) ; 0 lef_crs{iG}{iR+1,iC+1}(2)]; end      
                
                if ~sig && pvl>0.2 && ~any(strcmpi([cmb_nme{iR} cmb_nme{iC}],crcrcrc{iG})) &&  ~any(strcmpi([cmb_nme{iC} cmb_nme{iR}],crcrcrc{iG}))
                    rrr{iG}{end+1} = cmb_nme{iR};
                    ccc{iG}{end+1} = cmb_nme{iC};
                    crcrcrc{iG}{end+1} = [cmb_nme{iR} cmb_nme{iC}];
                end
                
            else
                lef_crs{iG}{iR+1,iC+1} = '-';
            end
            
        end
        
        crs_hld         = crosstab(yyy,zzz);
        if size(crs_hld,2)==1; crs_hld = [0 crs_hld(1) ; 0 crs_hld(2)]; end
        lef_bar{iG}(iR) = crs_hld(1,1)/crs_hld(1,2);
        
    end
    
    figure(); aaa{iG} = graph(rrr{iG},ccc{iG}); plot(aaa{iG});
    
end

figure()
bar(lef_bar{1})
title('Letter')
set(gca,'XTickLabel',cmb_nme)

figure()
bar(lef_bar{2})
title('Word')
set(gca,'XTickLabel',cmb_nme)

figure()
bar(lef_bar{3})
title('False-Font')
set(gca,'XTickLabel',cmb_nme)

%% Within Region Word Comparisons
% LHS %%%%%%%%%%%%%%%%%%%%%%%%%%%5
nme = fieldnames(reg);
cmb_nme = fieldnames(cmb_reg);

tst_nme = 'pap_wrd_600';
hms     = { 'lhs' };
col     = { [2 3] [12 13] [22 23] };
col_use = {[3 2] [1 2]};

tst.lhs = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/' tst_nme '/total/total_' tst_nme '_lhs_table_plot']);
tst.lhs = tst.lhs(2:end-1,:);
tmp_hld = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/' 'pap_con_600' '/total/total_' 'pap_con_600' '_lhs_table_plot']);
tst.lhs = [tst.lhs tmp_hld(2:end-1,[2 3])];

tst.rhs = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/' tst_nme '/total/total_' tst_nme '_rhs_table_plot']);
tst.rhs = tst.rhs(2:end-1,:);
tmp_hld = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/' 'pap_con_600' '/total/total_' 'pap_con_600' '_rhs_table_plot']);
tst.rhs = [tst.rhs tmp_hld(2:end-1,[2 3])];

% Pooled
wrd_hld_sig = cell(1,1);
wrd_hld_pvl = cell(1,1);
wrd_hld_crs = cell(1,1);
for iG = 1:numel(col_use)
    for iR = 1:numel(nme)
        
        reg_1st = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], reg.(nme{iR}),'uni',0)));
        reg_2nd = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], reg.(nme{iR}),'uni',0)));
        
        ele_num_1st = 0; tot_num_1st = 0; ele_num_2nd = 0; tot_num_2nd = 0;
        for iN = 1:numel(reg_1st)
            ele_num_1st = ele_num_1st + tst.(hms{1}){reg_1st(iN),col{col_use{iG}(1)}(1)};
            tot_num_1st = tot_num_1st + tst.(hms{1}){reg_1st(iN),col{col_use{iG}(1)}(2)};
            
            ele_num_2nd = ele_num_2nd + tst.(hms{1}){reg_2nd(iN),col{col_use{iG}(2)}(1)};
            tot_num_2nd = tot_num_2nd + tst.(hms{1}){reg_2nd(iN),col{col_use{iG}(2)}(2)};
        end
        
        yyy = [ ones(1,ele_num_1st)   ones(1,tot_num_1st-ele_num_1st)   ...
                ones(1,ele_num_2nd)*2 ones(1,tot_num_2nd-ele_num_2nd)*2 ];
        
        zzz = [ ones(1,ele_num_1st) ones(1,tot_num_1st-ele_num_1st)*2 ...
                ones(1,ele_num_2nd) ones(1,tot_num_2nd-ele_num_2nd)*2 ];
        
        [sig,pvl] = FisherExactTest(yyy,zzz);
        
        wrd_hld_sig{iR,1} = nme{iR};
        wrd_hld_pvl{iR,1} = nme{iR};
        wrd_hld_crs{iR,1} = nme{iR};
        
        wrd_hld_sig{iR,iG+1} = sig;
        wrd_hld_pvl{iR,iG+1} = pvl;        
        wrd_hld_crs{iR,iG+1} = crosstab(yyy,zzz);
        
        if size(wrd_hld_crs{iR,iG+1},2)==1; wrd_hld_crs{iR,iG+1} = [0 wrd_hld_crs{iR,iG+1}(1) ; 0 wrd_hld_crs{iR,iG+1}(2)]; end
        
    end
end

for iG = 1:numel(col)
    for iR = 1:numel(nme)
        bar_pct{iR,1}    = nme{iR};
        
        reg_1st = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], reg.(nme{iR}),'uni',0)));
        ele_num_1st = 0; tot_num_1st = 0;
        for iN = 1:numel(reg_1st)
            ele_num_1st = ele_num_1st + tst.(hms{1}){reg_1st(iN),col{iG}(1)};
            tot_num_1st = tot_num_1st + tst.(hms{1}){reg_1st(iN),col{iG}(2)};
        end
        
        bar_pct{iR,iG+1} = ele_num_1st/tot_num_1st;
    end
end

figure()
ttt = bar(cell2mat(bar_pct(:,2:4)));
ttt(1).FaceColor = rgb('purple'); ttt(2).FaceColor = rgb('red'); ttt(3).FaceColor = rgb('reddish grey');
set(gca,'XTickLabel',nme)

% Individual regions
wrd_hld_sig = cell(1,1);
wrd_hld_pvl = cell(1,1);
wrd_hld_crs = cell(1,1);
for iG = 1:numel(col_use)
    for iR = 1:numel(cmb_nme)
        
        reg_1st = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], cmb_reg.(cmb_nme{iR}),'uni',0)));
        reg_2nd = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], cmb_reg.(cmb_nme{iR}),'uni',0)));
        
        ele_num_1st = 0; tot_num_1st = 0; ele_num_2nd = 0; tot_num_2nd = 0;
        for iN = 1:numel(reg_1st)
            ele_num_1st = ele_num_1st + tst.(hms{1}){reg_1st(iN),col{col_use{iG}(1)}(1)};
            tot_num_1st = tot_num_1st + tst.(hms{1}){reg_1st(iN),col{col_use{iG}(1)}(2)};
            
            ele_num_2nd = ele_num_2nd + tst.(hms{1}){reg_2nd(iN),col{col_use{iG}(2)}(1)};
            tot_num_2nd = tot_num_2nd + tst.(hms{1}){reg_2nd(iN),col{col_use{iG}(2)}(2)};
        end
        
        yyy = [ ones(1,ele_num_1st)   ones(1,tot_num_1st-ele_num_1st)   ...
                ones(1,ele_num_2nd)*2 ones(1,tot_num_2nd-ele_num_2nd)*2 ];
        
        zzz = [ ones(1,ele_num_1st) ones(1,tot_num_1st-ele_num_1st)*2 ...
                ones(1,ele_num_2nd) ones(1,tot_num_2nd-ele_num_2nd)*2 ];
        
        [sig,pvl] = FisherExactTest(yyy,zzz);
        
        wrd_hld_sig{iR,1} = cmb_nme{iR};
        wrd_hld_pvl{iR,1} = cmb_nme{iR};
        wrd_hld_crs{iR,1} = cmb_nme{iR};
        
        wrd_hld_sig{iR,iG+1} = sig;
        wrd_hld_pvl{iR,iG+1} = pvl;        
        wrd_hld_crs{iR,iG+1} = crosstab(yyy,zzz);
        
        if size(wrd_hld_crs{iR,iG+1},2)==1; wrd_hld_crs{iR,iG+1} = [0 wrd_hld_crs{iR,iG+1}(1) ; 0 wrd_hld_crs{iR,iG+1}(2)]; end
        
    end
end

for iG = 1:numel(col)
    for iR = 1:numel(cmb_nme)
        bar_pct{iR,1}    = cmb_nme{iR};
        
        reg_1st = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], cmb_reg.(cmb_nme{iR}),'uni',0)));
        ele_num_1st = 0; tot_num_1st = 0;
        for iN = 1:numel(reg_1st)
            ele_num_1st = ele_num_1st + tst.(hms{1}){reg_1st(iN),col{iG}(1)};
            tot_num_1st = tot_num_1st + tst.(hms{1}){reg_1st(iN),col{iG}(2)};
        end
        
        bar_pct{iR,iG+1} = ele_num_1st/tot_num_1st;
    end
end

figure()
ttt = bar(cell2mat(bar_pct(:,2:4)));
ttt(1).FaceColor = rgb('purple'); ttt(2).FaceColor = rgb('red'); ttt(3).FaceColor = rgb('reddish grey');
set(gca,'XTickLabel',cmb_nme)

tot_wrd_hld_pvl = wrd_hld_pvl;
for iR = 1:numel(cmb_nme)
    for iC = 2:3
        if tot_wrd_hld_pvl{iR,iC} < 0.05
            tot_wrd_hld_pvl{iR,iC} = 9;
        elseif tot_wrd_hld_pvl{iR,iC} < 0.10
            tot_wrd_hld_pvl{iR,iC} = 3;
        elseif tot_wrd_hld_pvl{iR,iC} < 0.20
             tot_wrd_hld_pvl{iR,iC} = 1;
        else
             tot_wrd_hld_pvl{iR,iC} = 0;
        end
    end
end

% RHS %%%%%%%%%%%%%%%%%%%%%%%%%%%5
nme = fieldnames(reg);
cmb_nme = fieldnames(cmb_reg);

tst_nme = 'pap_wrd_600';
hms     = { 'rhs' };
col     = { [2 3] [12 13] [22 23] };
col_use = {[3 2] [1 2]};

tst.lhs = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/' tst_nme '/total/total_' tst_nme '_lhs_table_plot']);
tst.lhs = tst.lhs(2:end-1,:);
tmp_hld = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/' 'pap_con_600' '/total/total_' 'pap_con_600' '_lhs_table_plot']);
tst.lhs = [tst.lhs tmp_hld(2:end-1,[2 3])];

tst.rhs = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/' tst_nme '/total/total_' tst_nme '_rhs_table_plot']);
tst.rhs = tst.rhs(2:end-1,:);
tmp_hld = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/' 'pap_con_600' '/total/total_' 'pap_con_600' '_rhs_table_plot']);
tst.rhs = [tst.rhs tmp_hld(2:end-1,[2 3])];

% Pooled
wrd_hld_sig = cell(1,1);
wrd_hld_pvl = cell(1,1);
wrd_hld_crs = cell(1,1);
for iG = 1:numel(col_use)
    for iR = 1:numel(nme)
        
        reg_1st = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], reg.(nme{iR}),'uni',0)));
        reg_2nd = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], reg.(nme{iR}),'uni',0)));
        
        ele_num_1st = 0; tot_num_1st = 0; ele_num_2nd = 0; tot_num_2nd = 0;
        for iN = 1:numel(reg_1st)
            ele_num_1st = ele_num_1st + tst.(hms{1}){reg_1st(iN),col{col_use{iG}(1)}(1)};
            tot_num_1st = tot_num_1st + tst.(hms{1}){reg_1st(iN),col{col_use{iG}(1)}(2)};
            
            ele_num_2nd = ele_num_2nd + tst.(hms{1}){reg_2nd(iN),col{col_use{iG}(2)}(1)};
            tot_num_2nd = tot_num_2nd + tst.(hms{1}){reg_2nd(iN),col{col_use{iG}(2)}(2)};
        end
        
        yyy = [ ones(1,ele_num_1st)   ones(1,tot_num_1st-ele_num_1st)   ...
                ones(1,ele_num_2nd)*2 ones(1,tot_num_2nd-ele_num_2nd)*2 ];
        
        zzz = [ ones(1,ele_num_1st) ones(1,tot_num_1st-ele_num_1st)*2 ...
                ones(1,ele_num_2nd) ones(1,tot_num_2nd-ele_num_2nd)*2 ];
        
        [sig,pvl] = FisherExactTest(yyy,zzz);
        
        wrd_hld_sig{iR,1} = nme{iR};
        wrd_hld_pvl{iR,1} = nme{iR};
        wrd_hld_crs{iR,1} = nme{iR};
        
        wrd_hld_sig{iR,iG+1} = sig;
        wrd_hld_pvl{iR,iG+1} = pvl;        
        wrd_hld_crs{iR,iG+1} = crosstab(yyy,zzz);
        
        if size(wrd_hld_crs{iR,iG+1},2)==1; wrd_hld_crs{iR,iG+1} = [0 wrd_hld_crs{iR,iG+1}(1) ; 0 wrd_hld_crs{iR,iG+1}(2)]; end
        
    end
end

for iG = 1:numel(col)
    for iR = 1:numel(nme)
        bar_pct{iR,1}    = nme{iR};
        
        reg_1st = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], reg.(nme{iR}),'uni',0)));
        ele_num_1st = 0; tot_num_1st = 0;
        for iN = 1:numel(reg_1st)
            ele_num_1st = ele_num_1st + tst.(hms{1}){reg_1st(iN),col{iG}(1)};
            tot_num_1st = tot_num_1st + tst.(hms{1}){reg_1st(iN),col{iG}(2)};
        end
        
        bar_pct{iR,iG+1} = ele_num_1st/tot_num_1st;
    end
end

figure()
ttt = bar(cell2mat(bar_pct(:,2:4)));
ttt(1).FaceColor = rgb('purple'); ttt(2).FaceColor = rgb('red'); ttt(3).FaceColor = rgb('reddish grey');
set(gca,'XTickLabel',nme)

% Individual regions
wrd_hld_sig = cell(1,1);
wrd_hld_pvl = cell(1,1);
wrd_hld_crs = cell(1,1);
for iG = 1:numel(col_use)
    for iR = 1:numel(cmb_nme)
        
        reg_1st = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], cmb_reg.(cmb_nme{iR}),'uni',0)));
        reg_2nd = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], cmb_reg.(cmb_nme{iR}),'uni',0)));
        
        ele_num_1st = 0; tot_num_1st = 0; ele_num_2nd = 0; tot_num_2nd = 0;
        for iN = 1:numel(reg_1st)
            ele_num_1st = ele_num_1st + tst.(hms{1}){reg_1st(iN),col{col_use{iG}(1)}(1)};
            tot_num_1st = tot_num_1st + tst.(hms{1}){reg_1st(iN),col{col_use{iG}(1)}(2)};
            
            ele_num_2nd = ele_num_2nd + tst.(hms{1}){reg_2nd(iN),col{col_use{iG}(2)}(1)};
            tot_num_2nd = tot_num_2nd + tst.(hms{1}){reg_2nd(iN),col{col_use{iG}(2)}(2)};
        end
        
        yyy = [ ones(1,ele_num_1st)   ones(1,tot_num_1st-ele_num_1st)   ...
                ones(1,ele_num_2nd)*2 ones(1,tot_num_2nd-ele_num_2nd)*2 ];
        
        zzz = [ ones(1,ele_num_1st) ones(1,tot_num_1st-ele_num_1st)*2 ...
                ones(1,ele_num_2nd) ones(1,tot_num_2nd-ele_num_2nd)*2 ];
        
        [sig,pvl] = FisherExactTest(yyy,zzz);
        
        wrd_hld_sig{iR,1} = cmb_nme{iR};
        wrd_hld_pvl{iR,1} = cmb_nme{iR};
        wrd_hld_crs{iR,1} = cmb_nme{iR};
        
        wrd_hld_sig{iR,iG+1} = sig;
        wrd_hld_pvl{iR,iG+1} = pvl;        
        wrd_hld_crs{iR,iG+1} = crosstab(yyy,zzz);
        
        if size(wrd_hld_crs{iR,iG+1},2)==1; wrd_hld_crs{iR,iG+1} = [0 wrd_hld_crs{iR,iG+1}(1) ; 0 wrd_hld_crs{iR,iG+1}(2)]; end
        
    end
end

for iG = 1:numel(col)
    for iR = 1:numel(cmb_nme)
        bar_pct{iR,1}    = cmb_nme{iR};
        
        reg_1st = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], cmb_reg.(cmb_nme{iR}),'uni',0)));
        ele_num_1st = 0; tot_num_1st = 0;
        for iN = 1:numel(reg_1st)
            ele_num_1st = ele_num_1st + tst.(hms{1}){reg_1st(iN),col{iG}(1)};
            tot_num_1st = tot_num_1st + tst.(hms{1}){reg_1st(iN),col{iG}(2)};
        end
        
        bar_pct{iR,iG+1} = ele_num_1st/tot_num_1st;
    end
end

figure()
ttt = bar(cell2mat(bar_pct(:,2:4)));
ttt(1).FaceColor = rgb('purple'); ttt(2).FaceColor = rgb('red'); ttt(3).FaceColor = rgb('reddish grey');
set(gca,'XTickLabel',cmb_nme)

tot_wrd_hld_pvl = wrd_hld_pvl;
for iR = 1:numel(cmb_nme)
    for iC = 2:3
        if tot_wrd_hld_pvl{iR,iC} < 0.05
            tot_wrd_hld_pvl{iR,iC} = 9;
        elseif tot_wrd_hld_pvl{iR,iC} < 0.10
            tot_wrd_hld_pvl{iR,iC} = 3;
        elseif tot_wrd_hld_pvl{iR,iC} < 0.20
             tot_wrd_hld_pvl{iR,iC} = 1;
        else
             tot_wrd_hld_pvl{iR,iC} = 0;
        end
    end
end

%% Lex Overlap
frq_ovr = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/pap_lex_600/subjects/total' '/' 'pap_lex_600_plt']);
loc_ovr = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/electrode_location_files/total/output' '/' 'total_lhs_ecog']);
[~,frq_eff_lhs,~] = intersect(frq_ovr(:,2),loc_ovr(:,1));

lex_eff = cell2mat(frq_ovr(frq_eff_lhs,5)); % lex_eff = cell2mat(frq_ovr(2:end,5));
rep_eff = cell2mat(frq_ovr(frq_eff_lhs,3)); % rep_eff = cell2mat(frq_ovr(2:end,3));

sum(lex_eff & rep_eff) / sum(lex_eff);

%% Lex Locations
% LOAD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nme = fieldnames(reg);
cmb_nme = fieldnames(cmb_reg);

reg.tot_reg = cell(0);
for iN = 1:numel(nme)
    reg.tot_reg = [reg.tot_reg reg.(nme{iN})];
end
nme{end+1} = 'tot_reg';

cmb_reg.tot_reg = cell(0);
for iN = 1:numel(cmb_nme)
    cmb_reg.tot_reg = [cmb_reg.tot_reg cmb_reg.(cmb_nme{iN})];
end
cmb_nme{end+1} = 'tot_reg';

% COMPARING HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tst_nme = 'pap_lex_600';
hms     = {'lhs' 'rhs'};
col     = {[2 3] [12 13]};

tst.lhs = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/' tst_nme '/total/total_' tst_nme '_lhs_table_plot']);
tst.lhs = tst.lhs(2:end-1,:);

tst.rhs = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/' tst_nme '/total/total_' tst_nme '_rhs_table_plot']);
tst.rhs = tst.rhs(2:end-1,:);

cmb_hld = cell(1,1);
for iG = 1:numel(col)
    for iR = 1:numel(cmb_nme)
        
        reg_1st = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], cmb_reg.(cmb_nme{iR}),'uni',0)));
        reg_2nd = find(ismember(tst.(hms{2})(:,1),cellfun(@(x) [hms{2} '_' x], cmb_reg.(cmb_nme{iR}),'uni',0)));
        
        ele_num_1st = 0; tot_num_1st = 0; ele_num_2nd = 0; tot_num_2nd = 0;
        for iN = 1:numel(reg_1st)
            ele_num_1st = ele_num_1st + tst.(hms{1}){reg_1st(iN),col{iG}(1)};
            tot_num_1st = tot_num_1st + tst.(hms{1}){reg_1st(iN),col{iG}(2)};
            
            ele_num_2nd = ele_num_2nd + tst.(hms{2}){reg_2nd(iN),col{iG}(1)};
            tot_num_2nd = tot_num_2nd + tst.(hms{2}){reg_2nd(iN),col{iG}(2)};
        end
        
        yyy = [ ones(1,ele_num_1st)   ones(1,tot_num_1st-ele_num_1st)   ...
            ones(1,ele_num_2nd)*2 ones(1,tot_num_2nd-ele_num_2nd)*2 ];
        
        zzz = [ ones(1,ele_num_1st) ones(1,tot_num_1st-ele_num_1st)*2 ...
            ones(1,ele_num_2nd) ones(1,tot_num_2nd-ele_num_2nd)*2 ];
        
        [sig,pvl] = FisherExactTest(yyy,zzz);
        
        cmb_hld{iG}{iR,1} = cmb_nme{iR};
        cmb_hld{iG}{iR,2} = sig;
        cmb_hld{iG}{iR,3} = pvl;
        cmb_hld{iG}{iR,4} = crosstab(yyy,zzz);
        if size(cmb_hld{iG}{iR,4},2)==1; cmb_hld{iG}{iR,4} = [0 cmb_hld{iG}{iR,4}(1) ; 0 cmb_hld{iG}{iR,4}(2)]; end
        
    end
end

pct_hld = cat(3,cmb_hld{1}{:,4});
pct_bar = [squeeze(pct_hld(1,1,:))./(squeeze(pct_hld(1,2,:))+squeeze(pct_hld(1,1,:))) squeeze(pct_hld(2,1,:))./(squeeze(pct_hld(2,2,:))+squeeze(pct_hld(2,1,:)))];
figure()
bar(pct_bar)
set(gca,'XTickLabel',cmb_nme)
[strcat('lhs_',cmb_hld{1}(:,1)) num2cell(squeeze(pct_hld(1,1,:))) strcat('rhs_',cmb_hld{1}(:,1)) num2cell(squeeze(pct_hld(2,1,:)))]

pct_hld = cat(3,cmb_hld{2}{:,4});
pct_bar = [squeeze(pct_hld(1,1,:))./(squeeze(pct_hld(1,2,:))+squeeze(pct_hld(1,1,:))) squeeze(pct_hld(2,1,:))./(squeeze(pct_hld(2,2,:))+squeeze(pct_hld(2,1,:)))];
figure()
bar(pct_bar)
set(gca,'XTickLabel',cmb_nme)
[strcat('lhs_',cmb_hld{2}(:,1)) num2cell(squeeze(pct_hld(1,1,:))) strcat('rhs_',cmb_hld{2}(:,1)) num2cell(squeeze(pct_hld(2,1,:)))]

% COMPARING Within %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


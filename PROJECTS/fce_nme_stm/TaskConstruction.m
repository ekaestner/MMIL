clear; clc;

ovr_fld = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/iEEG/fce_nme_stm/Task';
img_fld = [ ovr_fld '/' 'cfd' '/' 'CFD Version 3.0' '/' 'Images'];

img_out_fld = [ ovr_fld '/' 'FN_stm_immediate_recall_v1' '/' 'stims' ];
lst_out_fld = [ ovr_fld '/' 'FN_stm_immediate_recall_v1' '/' 'lists' ]; 
ctr_out_fld = [ ovr_fld '/' 'FN_stm_immediate_recall_v1' '/' 'control' ]; 

lst_dly_out_fld = [ ovr_fld '/' 'FN_stm_delayed_recall_30m_v1' '/' 'lists' ]; 

num_lst     = 15;
num_lst_stm = 6;
rcl_enc_num = 3;
pix_max     = 800;

inp_tme = 3500; % i
rec_tme = 6000; % pl
brk_tme = 8000; % b

stm_col     = 1;
nme_col     = 2;
typ_col     = 3;
ord_col     = 4;
int_stm_col = 5;
blk_stm_col = 6;

lst_key     = { 'stimuli' 'name' 'trialtype' 'serialposition' 'interleave_stim' 'block_stim' };

%% Stimulation Decision
int_stm = { 'S' 'L' 'R' ; ...
            'L' 'R' 'S' ; ...
            'R' 'S' 'L' };

blk_stm = { 'S' 'L' 'R' ; ...
            'S' 'L' 'R' ; ...
            'S' 'L' 'R' };

int_stm = repmat( int_stm, num_lst_stm/size(int_stm,1),1);
blk_stm = repmat( blk_stm, num_lst_stm/size(blk_stm,1),1);
stm_ord = [ repmat(1:3,1,num_lst/size(blk_stm,2)) 1];
        
%% Load Faces
fce_hld = mmil_readtext([ ovr_fld '/' 'cfd' '/' 'CFD Version 3.0' '/' 'CFD_3.0_Norming_Data.csv']);
    fce_hld_col = fce_hld(1:6,:);
    fce_hld     = fce_hld(7:end,:);

% Split Ethnicity & Gender
fce_typ = { 'AF' 'AM' 'BF' 'BM' 'LF' 'LM' 'WF' 'WM' 'MF' 'MM' };
eth_typ = cellfun(@(x) {x(1)},fce_typ);

fce.col_nme = fce_hld_col;
for iFT = 1:numel(fce_typ)
    fce.(fce_typ{iFT}) = fce_hld(string_find(fce_hld(:,1),{[fce_typ{iFT} '-']}),:);
end

% Grab norming variables    
fce_nor_var = {};
fce_stm_var_col = unique( [string_find(fce.col_nme(2,:),{'Model'}) string_find(fce.col_nme(5,:),{'Model'})] );
    if numel(fce_stm_var_col)>1; error('naming column error'); end
fce_nor_var_col = { 'AgeRated' 'Attractive' 'Threatening' 'Unusual' 'LuminanceMedian' 'fWHR2' 'FaceShape' };
    fce_nor_var_col_num = numel(fce_nor_var_col);
fce_nor_var_col = unique( [string_find(fce.col_nme(2,:),fce_nor_var_col) string_find(fce.col_nme(5,:),fce_nor_var_col)] );
    if numel(fce_nor_var_col)>fce_nor_var_col_num; error('naming column error'); end

%% Load Names
nme_hld = mmil_readtext([ ovr_fld '/' 'Names_Popular_1990_nocomma.csv']);
    nme_hld_col = nme_hld(2,:);
    nme_hld     = nme_hld(3:end,:);
    
% Split Gender
nme.col_nme = [ nme_hld_col(2:3) 'Number' 'Gender'];
    nme.col_nme{2} = 'Frequency';
nme.M     = nme_hld(:,2:3);
    nme.M(:,end+1) = num2cell(cellfun(@numel,nme.M(:,1)));
    nme.M(:,end+1) = repmat({'M'},size(nme.M,1),1);
nme.F     = nme_hld(:,4:5);
    nme.F(:,end+1) = num2cell(cellfun(@numel,nme.F(:,1)));
    nme.F(:,end+1) = repmat({'F'},size(nme.F,1),1);
    
% Grab norming variables
nme_stm_var_col = string_find(nme.col_nme,{'Name'});
nme_nor_var_col = string_find(nme.col_nme,{ 'Gender' 'Number' 'Frequency' });

%% Make Lists
top_enc_lne = [ {'blank.png'} {'LEARN'}    {'b'} {0} {'S'} {'S'} repmat({''},1,numel(fce_nor_var_col)+numel(nme_nor_var_col))];
top_rcl_lne = [ {'blank.png'} {'REMEMBER'} {'b'} {0} {'S'} {'S'} repmat({''},1,numel(fce_nor_var_col)+numel(nme_nor_var_col))];

lst_key = [ lst_key nme.col_nme(nme_nor_var_col) fce.col_nme(2,fce_nor_var_col) ];

for iNL = 1:num_lst+1
    
    lst.(['L_' num2str(iNL)]) = cell((num_lst_stm+1)*(rcl_enc_num+1),numel(top_enc_lne));
    fce_img.(['L_' num2str(iNL)]) = cell(num_lst_stm,1);
    
    wht_ind = 0;
    while ~wht_ind
        eth_gen = eth_typ(randsample(numel(eth_typ),numel(eth_typ)));
        eth_gen = eth_gen(1:num_lst_stm);
        if sum(strcmpi(eth_gen,'W'))==2; wht_ind=1; end
    end
    
    lst_gen = [repmat({'M'},1,num_lst_stm/2) repmat({'F'},1,num_lst_stm/2)];
    lst_gen = lst_gen(randsample(num_lst_stm,num_lst_stm));
    
    fce_idn = cellfun( @(x,y) [x y],eth_gen,lst_gen,'uni',0);
    
    nme_use = cell(num_lst_stm,1);
    nme_use_nor = cell(num_lst_stm,numel(nme_nor_var_col));
    fce_use = cell(num_lst_stm,1);
    fce_use_nor = cell(num_lst_stm,numel(fce_nor_var_col));
    
    for iFS = 1:num_lst_stm
        
        % Choose Name        
        nme_use_ind        = randsample(1:size(nme.(lst_gen{iFS}),1),1);
        nme_use(iFS,1)     = upper(nme.(lst_gen{iFS})(nme_use_ind,nme_stm_var_col));
        nme_use_nor(iFS,:) = nme.(lst_gen{iFS})(nme_use_ind,nme_nor_var_col);     
        nme.(lst_gen{iFS})(nme_use_ind,:) = [];
        
        % Choose Face
        fce_use_ind        = randsample(1:size(fce.(fce_idn{iFS}),1),1);        
        fce_use(iFS,1)     = fce.(fce_idn{iFS})(fce_use_ind,fce_stm_var_col);  
        fce_use_nor(iFS,:) = fce.(fce_idn{iFS})(fce_use_ind,fce_nor_var_col);
        
        if strcmpi(fce_use{iFS,1}(1:3),'MM-') || strcmpi(fce_use{iFS,1}(1:3),'MF-')
            img_lod_fld = [ img_fld '/' 'CFD-MR'];
            pos_img     = dir(img_lod_fld);
            pos_img = {pos_img(:).name};
            pos_img = pos_img(string_find(pos_img,fce_use(iFS,1)));
            fce_use(iFS,1) = pos_img(string_find(pos_img,'-N.jpg'));
        else
            img_lod_fld = [ img_fld '/' 'CFD' '/' fce_use{iFS,1}];
            pos_img     = dir(img_lod_fld);
            fce_use(iFS,1) = {pos_img(string_find({pos_img(:).name},'-N.jpg')).name}; 
        end
        
        fce_img_fld.(['L_' num2str(iNL)]){iFS} = [img_lod_fld '/' fce_use{iFS,1}];
        fce_img_fle.(['L_' num2str(iNL)]){iFS} = fce_use{iFS,1};   
        
        fce.(fce_idn{iFS})(fce_use_ind,:) = [];
               
    end
    
     % Input list
     lst.(['L_' num2str(iNL)])(1,:)               = top_enc_lne;
     lst.(['L_' num2str(iNL)])(2:num_lst_stm+1,:) = [ fce_use nme_use repmat({'i'},num_lst_stm,1) num2cell(1:num_lst_stm)' int_stm(:,stm_ord(iNL)) blk_stm(:,stm_ord(iNL)) nme_use_nor fce_use_nor];
     
     % Recall/Encode lists
     for iRC = 1:rcl_enc_num
         lst_ind = ((num_lst_stm+1)*(iRC-1))+(num_lst_stm+2):((num_lst_stm+1)*(iRC-1))+(num_lst_stm*2+2);        
         lst.(['L_' num2str(iNL)])(lst_ind(1),:)               = top_rcl_lne;
         lst.(['L_' num2str(iNL)])(lst_ind(2:num_lst_stm+1),:) = [ fce_use nme_use repmat({'pl'},num_lst_stm,1) num2cell(1:num_lst_stm)' int_stm(:,stm_ord(iNL)) blk_stm(:,stm_ord(iNL)) nme_use_nor fce_use_nor];
         
         lst_ind_reo = lst_ind(2:num_lst_stm+1);
         
         % Flipping
         mve_rnd = randsample(3,1);
         
         if mve_rnd==1;     new_stm = int_stm(:,1);
         elseif mve_rnd==2; new_stm = int_stm(:,2);
         elseif mve_rnd==3; new_stm = int_stm(:,3);
         end 
         
         old_shm_ind = find(strcmpi(lst.(['L_' num2str(iNL)])(lst_ind_reo,5),'S')); new_shm_ind = find(strcmpi(new_stm,'S'));
         old_lft_ind = find(strcmpi(lst.(['L_' num2str(iNL)])(lst_ind_reo,5),'L')); new_lft_ind = find(strcmpi(new_stm,'L'));
         old_rgh_ind = find(strcmpi(lst.(['L_' num2str(iNL)])(lst_ind_reo,5),'R')); new_rgh_ind = find(strcmpi(new_stm,'R'));
         
         new_lst_ind_reo(new_shm_ind) = lst_ind_reo(old_shm_ind);
         new_lst_ind_reo(new_lft_ind) = lst_ind_reo(old_lft_ind);
         new_lst_ind_reo(new_rgh_ind) = lst_ind_reo(old_rgh_ind);

         % Implement flip
         lst.(['L_' num2str(iNL)])(lst_ind(2:num_lst_stm+1),:) = lst.(['L_' num2str(iNL)])(new_lst_ind_reo,:);     
         
         clear new_lst_ind_reo
     end
    
     if iNL==num_lst+1
         prc_lst = lst.(['L_' num2str(iNL)]);         
         lst = rmfield(lst,['L_' num2str(iNL)]);
         
         prc_lst(:,int_stm_col) = prc_lst(:,blk_stm_col);
         
         ord_hld = cell2mat(prc_lst(:,ord_col));
                 
         cell2csv([ lst_out_fld '/' 'practice_one.csv' ],prc_lst(ismember(ord_hld,[0 1 2 3]),:))
         cell2csv([ lst_out_fld '/' 'practice_two.csv' ],prc_lst(ismember(ord_hld,[0 4 5 6]),:))
         
     else
         cell2csv([ lst_out_fld '/' 'l' num2str(iNL) '-d1_stm.csv' ],lst.(['L_' num2str(iNL)]))
     end
end

nme_typ = {'M' 'F'};
nme_out = [];
for iNT = 1:numel(nme_typ)
    nme_out = [ nme_out ; nme.(nme_typ{iNT}) ];
end
cell2csv([ lst_out_fld '/' 'RemainingNames.csv' ],nme_out)

fce_out = [];
for iFT = 1:numel(fce_typ)
    fce_out = [ fce_out ; fce.(fce_typ{iFT}) ];
end
cell2csv([ lst_out_fld '/' 'RemainingFace.csv' ],fce_out)

cell2csv([ lst_out_fld '/' 'ListKey.csv' ],lst_key)

%% Caculate Distance between items
dst_hld     = lst;
dst_fld_nme = fieldnames(dst_hld);

for iD = 1:numel(dst_fld_nme)
    
   dst_out.(dst_fld_nme{iD}) = cell(size(dst_hld.(dst_fld_nme{iD}),1),3);
    
   for iST = 1:size(dst_hld.(dst_fld_nme{iD}),1)
       
        dst_out.(dst_fld_nme{iD})(iST,2) = dst_hld.(dst_fld_nme{iD})(iST,int_stm_col);
        dst_out.(dst_fld_nme{iD})(iST,3) = dst_hld.(dst_fld_nme{iD})(iST,typ_col);
        
        if strcmpi(dst_hld.(dst_fld_nme{iD})(iST,typ_col),'pl')
            
            nme_pst_ind = find(strcmpi(dst_hld.(dst_fld_nme{iD})(:,nme_col),dst_hld.(dst_fld_nme{iD}){iST,nme_col}));
            nme_pst_ind = nme_pst_ind(find((nme_pst_ind-iST)==0)-1);
            
            fce_pst_ind = find(strcmpi(dst_hld.(dst_fld_nme{iD})(:,stm_col),dst_hld.(dst_fld_nme{iD}){iST,stm_col}));
            fce_pst_ind = fce_pst_ind(find((fce_pst_ind-iST)==0)-1);            
            
            if nme_pst_ind ~= fce_pst_ind; error('Face/Name mismatch'); end
            
            typ_hld = dst_hld.(dst_fld_nme{iD})(nme_pst_ind+1:iST-1,typ_col);
            
            dst_out.(dst_fld_nme{iD}){iST,1} = sum(strcmpi(typ_hld,'i'))*inp_tme + ...
                                              sum(strcmpi(typ_hld,'pl'))*rec_tme + ...
                                              sum(strcmpi(typ_hld,'b'))*brk_tme;
            
        else
            dst_out.(dst_fld_nme{iD}){iST,1} = NaN;
        end        
   end
   
   lst.(['L_' num2str(iD)])(:,end+1) = dst_out.(dst_fld_nme{iD})(:,1);
   
end

lst_key = [ lst_key 'Distance_ms' ];

%% Check list control
% Distance between items %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
% Plot
plt_ind = 3:3:num_lst;
figure('Visible','off')
sub_plt_sqr = ceil(sqrt((numel(plt_ind))));

for iP = 1:numel(plt_ind)
    
    dta_use_hld = [];
    for iGR = 1:plt_ind(iP)
        dta_use_hld = [ dta_use_hld ; dst_out.(dst_fld_nme{iGR}) ];
    end
       
    sub_plt_ind = subplot(sub_plt_sqr,sub_plt_sqr,iP);
    
    fcfg = [];
    
    fcfg.xdt     = { 1 2 3 };
    fcfg.ydt     = { dta_use_hld(strcmpi(dta_use_hld(:,2),'S'),1)./1000 dta_use_hld(strcmpi(dta_use_hld(:,2),'L'),1)./1000 dta_use_hld(strcmpi(dta_use_hld(:,2),'R'),1)./1000 };
        fcfg.ydt = cellfun(@cell2mat,fcfg.ydt,'uni',0);

    pvl_hld = anova1( [ fcfg.ydt{1}(~isnan(fcfg.ydt{1})) fcfg.ydt{2}(~isnan(fcfg.ydt{2})) fcfg.ydt{3}(~isnan(fcfg.ydt{3})) ],[],'off');
        pvl_hld = num2str(roundsd(pvl_hld,2));
        pvl_hld = pvl_hld(2:end);
        
    fcfg.fce_col = { rgb('dark grey') rgb('dark teal') rgb('dark mauve')};
    fcfg.edg_col = repmat({rgb('black')},1,numel(fcfg.xdt));
    
    fcfg.xlb = { 'Sham' 'Left' 'Right' };
    fcfg.ylb = { 'ISI (s)'  };
        
    fcfg.xlm = [ 0 4 ];
    
    fcfg.sbp = sub_plt_ind;
    fcfg.ttl = [ num2str(plt_ind(iP)) 'G p=' pvl_hld];
    
    ejk_scatter(fcfg)
        
end
    
print(gcf,[ ctr_out_fld '/' 'ISI_Distance.png'],'-dpng')
close all

end

% nuisance comparison across groups  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
    
    for iNV = 8:numel(lst_key)-1
                
        % Plot
        plt_ind = 3:3:num_lst;
        figure('Visible','off')
        sub_plt_sqr = ceil(sqrt((numel(plt_ind))));
        
        for iP = 1:numel(plt_ind)
            
            dta_use_hld = [];
            for iGR = 1:plt_ind(iP)
                dta_use_hld = [ dta_use_hld ; lst.(['L_' num2str(iGR)])(strcmpi(lst.(['L_' num2str(iGR)])(:,typ_col),'i'),[ iNV int_stm_col]) ];
            end
            
            sub_plt_ind = subplot(sub_plt_sqr,sub_plt_sqr,iP);
            
            fcfg = [];
            
            fcfg.xdt     = { 1 2 3 };
            fcfg.ydt     = { dta_use_hld(strcmpi(dta_use_hld(:,2),'S'),1) dta_use_hld(strcmpi(dta_use_hld(:,2),'L'),1) dta_use_hld(strcmpi(dta_use_hld(:,2),'R'),1) };
            fcfg.ydt = cellfun(@cell2mat,fcfg.ydt,'uni',0);
            
            pvl_hld = anova1( [ fcfg.ydt{1}(~isnan(fcfg.ydt{1})) fcfg.ydt{2}(~isnan(fcfg.ydt{2})) fcfg.ydt{3}(~isnan(fcfg.ydt{3})) ],[],'off');
            pvl_hld = num2str(roundsd(pvl_hld,2));
            pvl_hld = pvl_hld(2:end);
            
            fcfg.fce_col = { rgb('dark grey') rgb('dark teal') rgb('dark mauve')};
            fcfg.edg_col = repmat({rgb('black')},1,numel(fcfg.xdt));
            
            fcfg.xlb = { 'Sham' 'Left' 'Right' };
            fcfg.ylb = { lst_key{iNV}  };
            
            fcfg.xlm = [ 0 4 ];
            
            fcfg.sbp = sub_plt_ind;
            fcfg.ttl = [ num2str(plt_ind(iP)) 'G p=' pvl_hld];
            
            ejk_scatter(fcfg)
            
        end
        
        print(gcf,[ ctr_out_fld '/' lst_key{iNV} '.png'],'-dpng')
        close all
        
    end
end

% nuisance comparison across encoding serial position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nuisance comparison across test serial position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Copy over Faces
fce_fld_nme = fieldnames(fce_img_fld);
for iGR = 1:numel(fce_fld_nme)
        
    for iFS = 1:numel(fce_img_fld.(fce_fld_nme{iGR}))
        
        img_hld = imread(fce_img_fld.(fce_fld_nme{iGR}){iFS});
            img_rsz_prp = min(pix_max./size(img_hld,1:2));
            img_hld = imresize(img_hld,img_rsz_prp);
        imwrite(img_hld,[ img_out_fld '/' fce_img_fle.(fce_fld_nme{iGR}){iFS} ]);
        
    end
    
end

%% Modify lists for different stimulation paradigms
lst_fle_nme = dir(lst_out_fld);
    lst_fle_nme = {lst_fle_nme(:).name};
    prc_lst_nme = lst_fle_nme(string_find(lst_fle_nme,{'practice_'}));
    lst_fle_nme = lst_fle_nme(string_find(lst_fle_nme,{'_stm.csv'}));
    
for iL = 1:numel(lst_fle_nme)
    
    lst_hld = mmil_readtext( [ lst_out_fld '/' lst_fle_nme{iL} ]);
        lst_hld(cellfun( @isnumeric, lst_hld(:,typ_col)),typ_col) = {'i'};
        
    lst_hld = lst_hld(1:num_lst_stm+1,:);
    lst_hld{1,nme_col} = 'REMEMBER';
    lst_hld(strcmpi(lst_hld(:,typ_col),'i'),typ_col) = {'pl2'};
    lst_hld(:,int_stm_col) = repmat({'S'},size(lst_hld,1),1);
    
    cell2csv([ lst_dly_out_fld '/' strrep(lst_fle_nme{iL},'_stm.','_delayed_stm.') ],lst_hld)
        
end

%% Modify lists for different stimulation paradigms
lst_fle_nme = dir(lst_out_fld);
    lst_fle_nme = {lst_fle_nme(:).name};
    prc_lst_nme = lst_fle_nme(string_find(lst_fle_nme,{'practice_'}));
    lst_fle_nme = lst_fle_nme(string_find(lst_fle_nme,{'_stm.csv'}));

for iL = 1:numel(lst_fle_nme)
    
    lst_hld = mmil_readtext( [ lst_out_fld '/' lst_fle_nme{iL} ]);
        lst_hld(cellfun( @isnumeric, lst_hld(:,typ_col)),typ_col) = {'i'};
    
    lst_hld( strcmpi(lst_hld(:,typ_col),'i') & strcmpi(lst_hld(:,int_stm_col),'L'),int_stm_col)  = {'L_lrn'};
    lst_hld( strcmpi(lst_hld(:,typ_col),'i') & strcmpi(lst_hld(:,int_stm_col),'R'),int_stm_col)  = {'R_lrn'};
    lst_hld( strcmpi(lst_hld(:,typ_col),'pl') & strcmpi(lst_hld(:,int_stm_col),'L'),int_stm_col) = {'L_imt'};
    lst_hld( strcmpi(lst_hld(:,typ_col),'pl') & strcmpi(lst_hld(:,int_stm_col),'R'),int_stm_col) = {'R_imt'};
    lst_hld( string_find(lst_hld(:,int_stm_col),'S'),int_stm_col) = {'S'};
    
    lst_hld( strcmpi(lst_hld(:,typ_col),'i') & strcmpi(lst_hld(:,blk_stm_col),'L'),blk_stm_col)  = {'L_lrn'};
    lst_hld( strcmpi(lst_hld(:,typ_col),'i') & strcmpi(lst_hld(:,blk_stm_col),'R'),blk_stm_col)  = {'R_lrn'};
    lst_hld( strcmpi(lst_hld(:,typ_col),'pl') & strcmpi(lst_hld(:,blk_stm_col),'L'),blk_stm_col) = {'L_imt'};
    lst_hld( strcmpi(lst_hld(:,typ_col),'pl') & strcmpi(lst_hld(:,blk_stm_col),'R'),blk_stm_col) = {'R_imt'};
    lst_hld( string_find(lst_hld(:,blk_stm_col),'S'),blk_stm_col) = {'S'};    
    
    cell2csv([ lst_out_fld '/' lst_fle_nme{iL} ],lst_hld)
    
end


for iL = 1:numel(prc_lst_nme)
    
    lst_hld = mmil_readtext( [ lst_out_fld '/' prc_lst_nme{iL} ]);
        lst_hld(cellfun( @isnumeric, lst_hld(:,typ_col)),typ_col) = {'i'};
end

%% Make Mega list
lst_fle_nme = dir(lst_out_fld);
    lst_fle_nme = {lst_fle_nme(:).name};
    lst_fle_nme = lst_fle_nme(string_find(lst_fle_nme,{'_stm.csv'}));

meg_lst = [];
for iL = 1:numel(lst_fle_nme)
    
    lst_hld = mmil_readtext( [ lst_out_fld '/' lst_fle_nme{iL} ]);
        lst_hld(cellfun( @isnumeric, lst_hld(:,typ_col)),typ_col) = {'i'};
      
    brk_add = lst_hld(1,:);
        brk_add{1,2} = 'REST BREAK';
        
    meg_lst = [ meg_lst ; lst_hld ; brk_add ; brk_add];
        
end
cell2csv([ lst_out_fld '/' 'mega_list.csv' ],meg_lst);



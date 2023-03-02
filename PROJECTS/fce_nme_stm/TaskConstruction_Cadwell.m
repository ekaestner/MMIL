clear; clc;

stm_fld = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Methods/DataCollection/Stimulation/Task/Stimuli/';
fce_fld = [ stm_fld '/' 'cfd' '/' 'CFD Version 3.0'];


ovr_fld = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Methods/DataCollection/Stimulation/Task_v2/'; 
enc_fld = [ ovr_fld '/' 'FN_stm_v1' '/' ];    ejk_chk_dir(enc_fld);
ret_fld = [ ovr_fld '/' 'FN_stm_delayed_recall_v1' ];  ejk_chk_dir(ret_fld);

enc_lst_out_fld = [ enc_fld '/' 'lists' ];   ejk_chk_dir(enc_lst_out_fld);
enc_img_out_fld = [ enc_fld '/' 'stims' ];   ejk_chk_dir(enc_img_out_fld);
enc_ctr_out_fld = [ enc_fld '/' 'control' ]; ejk_chk_dir(enc_ctr_out_fld);
enc_scr_out_fld = [ enc_fld '/' 'scoring' ]; ejk_chk_dir(enc_scr_out_fld);

ret_out_fld     = [ enc_fld '/' 'lists' ];   ejk_chk_dir(ret_out_fld);
ret_img_out_fld = [ enc_fld '/' 'stims' ];   ejk_chk_dir(ret_img_out_fld);

stm_col     = 1; % 
nme_col     = 2; % 
typ_col     = 3; % 
ord_col     = 4; % 
int_stm_col = 5; % 
blk_stm_col = 6; % 

lst_key     = { 'stimuli' 'name' 'trialtype' 'serialposition' 'interleave_stim' 'block_stim' };

%% Stimuli setup
num_lst     = 15; % Number of lists to generate
num_lst_stm = 5; % Number of target items per slist
rcl_enc_num = 3; % Number of recall/learn trials following initial encode

enc_inp_tme = 3500; % i: encode trials
enc_rec_tme = 6000; % pl: encode recall/learn trials
enc_brk_tme = 8000; % b: encode break trials

ret_rcl_tme = 6000; % dpl: retrieval recall trials
ret_rec_tme = 6000; % dpl: retrieval recognition trials

pix_max     = 800; % Size of stimuli

%% Stimulation setup
int_stm = { 'L' 'R' 'L' ; ...
            'R' 'L' 'R' ; ...
            'S' 'R' 'L' ; ...
            'L' 'S' 'R' ; ...
            'R' 'L' 'S' };


blk_stm = { 'S' 'L' 'R' };

int_stm = repmat( int_stm, num_lst_stm/size(int_stm,1),1);
blk_stm = repmat( blk_stm, num_lst_stm/size(blk_stm,1),1);
stm_ord = [ repmat(1:3,1,num_lst/size(blk_stm,2)) 1 2];
        
%% Load Faces
fce_hld = mmil_readtext([ fce_fld '/' 'CFD_3.0_Norming_Data.csv']);
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
nme_hld = mmil_readtext([ stm_fld '/' 'Names_Popular_1990_nocomma.csv']);
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

%% ENCODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make Encoding Lists
% Instruction trials
top_enc_lne = [ {'blank.png'} {'LEARN'}    {'b'} {0} {'S'} {'S'} repmat({''},1,numel(fce_nor_var_col)+numel(nme_nor_var_col))];
top_rcl_lne = [ {'blank.png'} {'REMEMBER'} {'b'} {0} {'S'} {'S'} repmat({''},1,numel(fce_nor_var_col)+numel(nme_nor_var_col))];

% Add on column labels
lst_key = [ lst_key nme.col_nme(nme_nor_var_col) fce.col_nme(2,fce_nor_var_col) ];

% Loop through Encoding List Creation
for iNL = 1:num_lst+2
    
    % Initialize Lists & Stimuli
    lst.(['L_' num2str(iNL,'%02d')]) = cell((num_lst_stm+1)*(rcl_enc_num+1),numel(top_enc_lne));
    fce_img.(['L_enc_' num2str(iNL,'%02d')]) = cell(num_lst_stm,1);
    
    % Balance White selection
    wht_ind = 0;
    while ~wht_ind
        eth_gen = eth_typ(randsample(numel(eth_typ),numel(eth_typ)));
        eth_gen = eth_gen(1:num_lst_stm);
        if sum(strcmpi(eth_gen,'W'))==2; wht_ind=1; end
    end
    
    % Balance Gender selection
    gen_ind = 0;
    while ~gen_ind
        lst_gen = [repmat({'M'},1,num_lst_stm) repmat({'F'},1,num_lst_stm)];
        lst_gen = lst_gen(randsample(num_lst_stm*2,num_lst_stm));
        if abs(sum(strcmpi(lst_gen,'M'))-sum(strcmpi(lst_gen,'F')))<=1; gen_ind=1; end
    end
    
    % Combine Race & Gender
    fce_idn = cellfun( @(x,y) [x y],eth_gen,lst_gen,'uni',0);
    
    % Initialize Stimuli & Normative values
    nme_use = cell(num_lst_stm,1);
    nme_use_nor = cell(num_lst_stm,numel(nme_nor_var_col));
    fce_use = cell(num_lst_stm,1);
    fce_use_nor = cell(num_lst_stm,numel(fce_nor_var_col));
    
    % Loop through trials
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
        
        % Find Face File Name
        if strcmpi(fce_use{iFS,1}(1:3),'MM-') || strcmpi(fce_use{iFS,1}(1:3),'MF-')
            img_lod_fld = [ fce_fld '/' 'Images' '/' 'CFD-MR'];
            pos_img     = dir(img_lod_fld);
            pos_img = {pos_img(:).name};
            pos_img = pos_img(string_find(pos_img,fce_use(iFS,1)));
            fce_use(iFS,1) = pos_img(string_find(pos_img,'-N.jpg'));
        else
            img_lod_fld = [ fce_fld '/' 'Images' '/' 'CFD' '/' fce_use{iFS,1}];
            pos_img     = dir(img_lod_fld);
            fce_use(iFS,1) = {pos_img(string_find({pos_img(:).name},'-N.jpg')).name}; 
        end
        
        fce_img_fld.(['L_enc_' num2str(iNL,'%02d')]){iFS} = [img_lod_fld '/' fce_use{iFS,1}];
        fce_img_fle.(['L_enc_' num2str(iNL,'%02d')]){iFS} = fce_use{iFS,1};   
        
        fce.(fce_idn{iFS})(fce_use_ind,:) = []; % Remove face from circulation
               
    end
    
     % Input list
     lst.(['L_' num2str(iNL,'%02d')])(1,:)               = top_enc_lne;
     lst.(['L_' num2str(iNL,'%02d')])(2:num_lst_stm+1,:) = [ fce_use nme_use repmat({'i'},num_lst_stm,1) num2cell(1:num_lst_stm)' int_stm(:,stm_ord(iNL)) blk_stm(:,stm_ord(iNL)) nme_use_nor fce_use_nor];
     
     % Recall/Encode lists
     for iRC = 1:rcl_enc_num
         lst_ind = ((num_lst_stm+1)*(iRC-1))+(num_lst_stm+2):((num_lst_stm+1)*(iRC-1))+(num_lst_stm*2+2);        
         lst.(['L_' num2str(iNL,'%02d')])(lst_ind(1),:)               = top_rcl_lne;
         lst.(['L_' num2str(iNL,'%02d')])(lst_ind(2:num_lst_stm+1),:) = [ fce_use nme_use repmat({'pl'},num_lst_stm,1) num2cell(1:num_lst_stm)' int_stm(:,stm_ord(iNL)) blk_stm(:,stm_ord(iNL)) nme_use_nor fce_use_nor];
         
         lst_ind_reo = lst_ind(2:num_lst_stm+1);
         
         % Flipping based on interleaved stimulation order
         mve_rnd = randsample(3,1);
                  
         if mve_rnd==1;     new_stm = int_stm(:,1);
         elseif mve_rnd==2; new_stm = int_stm(:,2);
         elseif mve_rnd==3; new_stm = int_stm(:,3);
         end
         
         old_shm_ind = find(strcmpi(lst.(['L_' num2str(iNL,'%02d')])(lst_ind_reo,5),'S')); new_shm_ind = find(strcmpi(new_stm,'S'));
         old_lft_ind = find(strcmpi(lst.(['L_' num2str(iNL,'%02d')])(lst_ind_reo,5),'L')); new_lft_ind = find(strcmpi(new_stm,'L'));
         old_rgh_ind = find(strcmpi(lst.(['L_' num2str(iNL,'%02d')])(lst_ind_reo,5),'R')); new_rgh_ind = find(strcmpi(new_stm,'R'));
         
         new_lst_ind_reo(new_shm_ind) = lst_ind_reo(old_shm_ind);
         new_lst_ind_reo(new_lft_ind) = lst_ind_reo(old_lft_ind);
         new_lst_ind_reo(new_rgh_ind) = lst_ind_reo(old_rgh_ind);

         % Implement flip
         lst.(['L_' num2str(iNL,'%02d')])(lst_ind(2:num_lst_stm+1),:) = lst.(['L_' num2str(iNL,'%02d')])(new_lst_ind_reo,:);     
         
         clear new_lst_ind_reo
     end
    
     if iNL>num_lst
         % Create Practice Lists
         prc_lst = lst.(['L_' num2str(iNL,'%02d')]);         
                  
         prc_lst(:,int_stm_col) = prc_lst(:,blk_stm_col);
         
         ord_hld = cell2mat(prc_lst(:,ord_col));
                 
         cell2csv([ enc_lst_out_fld '/' 'practice_one.csv' ],prc_lst(ismember(ord_hld,[0 1 2 ]),:))
         cell2csv([ enc_lst_out_fld '/' 'practice_two.csv' ],prc_lst(ismember(ord_hld,[0 4 5 ]),:))
         
     else
         % Save List as-is
         cell2csv([ enc_lst_out_fld '/' 'l' num2str(iNL,'%02d') '-d1_stm.csv' ],lst.(['L_' num2str(iNL,'%02d')]))
     end
end

% Remove Practice Lists from Lists
lst = rmfield(lst,['L_' num2str(num_lst+1,'%02d')]);
lst = rmfield(lst,['L_' num2str(num_lst+2,'%02d')]);

% Save out the names that were not used
nme_typ = {'M' 'F'};
nme_out = [];
for iNT = 1:numel(nme_typ)
    nme_out = [ nme_out ; nme.(nme_typ{iNT}) ];
end
cell2csv([ enc_lst_out_fld '/' 'RemainingNames.csv' ],nme_out)

% Save out the Faces that were not used
fce_out = [];
for iFT = 1:numel(fce_typ)
    fce_out = [ fce_out ; fce.(fce_typ{iFT}) ];
end
cell2csv([ enc_lst_out_fld '/' 'RemainingFace.csv' ],fce_out)

% Save out column names
cell2csv([ enc_lst_out_fld '/' 'ListKey.csv' ],lst_key)

%% Make Delayed Recall/Recognition lists
% top lines
top_dly_rcl_lne = [ {'blank.png'} {'REMEMBER'}  {'b'} {0} {'S'} {'S'} repmat({''},1,numel(fce_nor_var_col)+numel(nme_nor_var_col))];
top_dly_rec_lne = [ {'blank.png'} {'RECOGNIZE'} {'b'} {0} {'S'} {'S'} repmat({''},1,numel(fce_nor_var_col)+numel(nme_nor_var_col))];

% Get list names
lst_nme = dir([ enc_lst_out_fld '/' '*.csv' ]);
lst_nme = {lst_nme(:).name};
lst_nme(string_find(lst_nme,'Remaining')) = [];
lst_nme(string_find(lst_nme,'ListKey')) = [];

% Get remaining names
rem_nme = mmil_readtext([ enc_lst_out_fld '/' 'RemainingNames.csv' ]);
for iNT = 1:numel(nme_typ)
    nme.(nme_typ{iNT}) = rem_nme(strcmpi(rem_nme(:,end),nme_typ{iNT}),:);
end

% Get remaining faces
rem_fce = mmil_readtext([ enc_lst_out_fld '/' 'RemainingFace.csv' ]);
for iFT = 1:numel(fce_typ)
    fce.(fce_typ{iFT}) = rem_fce(string_find(rem_fce(:,1),fce_typ{iFT}),:);
end

% Loop through lists
for iL = 1:numel(lst_nme)
    
    lst_hld = mmil_readtext([ enc_lst_out_fld '/' lst_nme{iL}]);
        
    % Find original list of names
    fst_lne_ind = find(strcmpi(lst_hld(:,nme_col),'LEARN'));
    scd_lne_ind = find(strcmpi(lst_hld(:,nme_col),'REMEMBER'),1);
    stm_hld     = lst_hld(fst_lne_ind+1:scd_lne_ind-1,:);
    
    % RECALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stm_hld(:,typ_col) = repmat({'delayed_recall'},size(stm_hld,1),1);

    % Save out
    cell2csv([ ret_out_fld '/' lst_nme{iL}(1:end-4) '_' 'recall' '.csv'],[ top_dly_rcl_lne ; stm_hld ])
    
    % RECOGNITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize Lists & Stimuli
    lst.(['L_' num2str(iNL,'%02d')]) = cell(size(stm_hld,1),numel(top_enc_lne));
    fce_img.(['L_rec_' num2str(iL,'%02d')]) = cell(size(stm_hld,1),1);
    
    % Get Race/Gender & Initialize new stimuli
    fce_idn = cellfun(@(x) x(2),regexp(stm_hld(:,stm_col),'-','split'))';
    
    % Initialize Stimuli & Normative values
    nme_use = cell(size(stm_hld,1),1);
    nme_use_nor = cell(size(stm_hld,1),numel(nme_nor_var_col));
    fce_use = cell(size(stm_hld,1),1);
    fce_use_nor = cell(size(stm_hld,1),numel(fce_nor_var_col));
    
    % Choose new Stimuli by looping through
    for iFS = 1:size(stm_hld,1)
        
        % Choose Name        
        nme_use_ind        = randsample(1:size(nme.(lst_gen{iFS}),1),1);
        nme_use(iFS,1)     = upper(nme.(lst_gen{iFS})(nme_use_ind,nme_stm_var_col));
        nme_use_nor(iFS,:) = nme.(lst_gen{iFS})(nme_use_ind,nme_nor_var_col);     
        nme.(lst_gen{iFS})(nme_use_ind,:) = [];
        
        % Choose Face
        fce_use_ind        = randsample(1:size(fce.(fce_idn{iFS}),1),1);        
        fce_use(iFS,1)     = fce.(fce_idn{iFS})(fce_use_ind,fce_stm_var_col);  
        fce_use_nor(iFS,:) = fce.(fce_idn{iFS})(fce_use_ind,fce_nor_var_col);
        
        % Find Face File Name
        if strcmpi(fce_use{iFS,1}(1:3),'MM-') || strcmpi(fce_use{iFS,1}(1:3),'MF-')
            img_lod_fld = [ fce_fld '/' 'Images' '/' 'CFD-MR'];
            pos_img     = dir(img_lod_fld);
            pos_img = {pos_img(:).name};
            pos_img = pos_img(string_find(pos_img,fce_use(iFS,1)));
            fce_use(iFS,1) = pos_img(string_find(pos_img,'-N.jpg'));
        else
            img_lod_fld = [ fce_fld '/' 'Images' '/' 'CFD' '/' fce_use{iFS,1}];
            pos_img     = dir(img_lod_fld);
            fce_use(iFS,1) = {pos_img(string_find({pos_img(:).name},'-N.jpg')).name}; 
        end
        
        fce_img_fld.(['L_rec_' num2str(iL,'%02d')]){iFS} = [img_lod_fld '/' fce_use{iFS,1}];
        fce_img_fle.(['L_rec_' num2str(iL,'%02d')]){iFS} = fce_use{iFS,1};   
        
        fce.(fce_idn{iFS})(fce_use_ind,:) = []; % Remove face from circulation
               
    end
    
    % Get additional stimuli
    lst.(['L_' num2str(iL,'%02d')]) = [ fce_use nme_use repmat({'i'},size(stm_hld,1),1) num2cell(1:size(stm_hld,1))' int_stm(1:size(stm_hld,1),stm_ord(iL)) blk_stm(1:size(stm_hld,1),stm_ord(iL)) nme_use_nor fce_use_nor];
     
    % Add in additional stimuli to lists
    new_old_lst = [ stm_hld                  repmat({'OLD'},size(stm_hld,1),1) ; ...
                    lst.(['L_' num2str(iL,'%02d')]) repmat({'NOVEL'},size(stm_hld,1),1) ];
    new_old_lst(:,typ_col) = repmat({'delayed_recognition'},size(new_old_lst,1),1);
            
    % Mix up & Check distances
    run_mix = 1;
    while run_mix
        lst_try = new_old_lst(randsample(size(new_old_lst,1),size(new_old_lst,1)),:);
        dff_nov = diff(find(strcmpi(lst_try(:,end),'NOVEL')));
        dff_old = diff(find(strcmpi(lst_try(:,end),'OLD')));
        if all(dff_nov<4) && all(dff_old<4) && ...
                ~all(dff_nov==1) && ~all(dff_old==1) && ...
                (sum(dff_nov==1) + sum(dff_old==1)) < 4 && ...
                (sum(dff_nov==3) + sum(dff_old==3)) < 2
            run_mix = 0;
            new_old_lst = lst_try;
        end
    end
       
    % Save out
    cell2csv([ ret_out_fld '/' lst_nme{iL}(1:end-4) '_' 'recognition' '.csv'],[ top_dly_rec_lne {''}; new_old_lst])
    
end

%% Make manual scoring tables
% Get list names
lst_nme = dir([ enc_lst_out_fld '/' '*.csv' ]);
lst_nme = {lst_nme(:).name};
lst_nme = lst_nme([ string_find(lst_nme,'practice_one.csv')  string_find(lst_nme,'practice_two.csv') string_find(lst_nme,'d1_stm.csv') ] );

% Initialize Table
tbl_out = cell( ((num_lst+2) * num_lst_stm) + ((num_lst+2)*2), (rcl_enc_num*3) + 2 );
row_ind = 2;
for iL = 1:numel(lst_nme)
    
    col_ind = 1;
    
    tbl_out{row_ind-1,col_ind} = lst_nme{iL}(1:end-4);
    
    % Find original list of names
    lst_hld = mmil_readtext([ enc_lst_out_fld '/' lst_nme{iL}]);
    scd_lne_ind = [ find(strcmpi(lst_hld(:,nme_col),'REMEMBER')) ; size(lst_hld,1)+1];
        
    % Blocks
    row_use = row_ind:row_ind+numel(2:scd_lne_ind(1)-1)-1;
    for iB = 1:rcl_enc_num        
        tbl_out(row_use,col_ind) = lst_hld(scd_lne_ind(iB)+1:scd_lne_ind(iB+1)-1,nme_col);
        col_ind = col_ind+3;
    end
    
    % Delayed Recall
    tbl_out(row_use,col_ind) = lst_hld(2:scd_lne_ind(1)-1,nme_col);
    
    % Iterate rows
    row_ind = row_ind + numel(2:scd_lne_ind(1)-1) + 2;
end

cell2csv([enc_scr_out_fld '/' 'test_key.csv' ],tbl_out);

%% Copy over Faces - NEED TO ADD RECOGNITION FACES
% Encoding
fce_fld_nme = fieldnames(fce_img_fld);
for iGR = 1:numel(fce_fld_nme)        
    for iFS = 1:numel(fce_img_fld.(fce_fld_nme{iGR}))        
        img_hld = imread(fce_img_fld.(fce_fld_nme{iGR}){iFS});
            img_rsz_prp = min(pix_max./size(img_hld,1:2));
            img_hld = imresize(img_hld,img_rsz_prp);
        imwrite(img_hld,[ enc_img_out_fld '/' fce_img_fle.(fce_fld_nme{iGR}){iFS} ]);
    end    
end

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
            
            dst_out.(dst_fld_nme{iD}){iST,1} = sum(strcmpi(typ_hld,'i'))*enc_inp_tme + ...
                                              sum(strcmpi(typ_hld,'pl'))*enc_rec_tme + ...
                                              sum(strcmpi(typ_hld,'b'))*enc_brk_tme;
            
        else
            dst_out.(dst_fld_nme{iD}){iST,1} = NaN;
        end        
   end
   
   lst.(['L_' num2str(iD,'%02d')])(:,end+1) = dst_out.(dst_fld_nme{iD})(:,1);
   
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
                dta_use_hld = [ dta_use_hld ; lst.(['L_' num2str(iGR,'%02d')])(strcmpi(lst.(['L_' num2str(iGR,'%02d')])(:,typ_col),'i'),[ iNV int_stm_col]) ];
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

% nuisance comparison across encoding serial position 

% nuisance comparison across test serial position 



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

%% Make recall lists



%% Make recognition lists









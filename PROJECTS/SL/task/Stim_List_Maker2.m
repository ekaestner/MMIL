clear; clc;

%% Make lists for SL
con = {'P' 'K' 'D' 'T' 'F' 'Z' 'S' 'W' 'L' 'R' 'V' 'Y' 'M' 'N' 'H'};
vow = {'AW' 'OO' 'UR' 'OW' 'AY' 'EH' 'EE' 'IH' 'UH'};

mes_con = [1 1 1 1 2 2 2 3 3 3 5 5 6 6 0];
mes_vow = [3 3 3 4 4 4 5 5 0];

vc_wrd = mmil_readtext('/home/ekaestne/Desktop/Tasks/VC_words/VC_wordlist_minusA.csv');
vc_wrd_cut = vc_wrd(1:end-2,1:2);

cv_wrd = mmil_readtext('/home/ekaestne/Desktop/Tasks/CV_words/CV_wordlist.csv');
cv_wrd_cut = cv_wrd(2:end,[1 3]);

wrd = {};
% Create Visual Stimuli

spk_num_mtc = repmat([1 2],1,numel(con)*numel(vow));
spk_num_mtc = spk_num_mtc(randperm(numel(spk_num_mtc)-1) + 1);

%% Visual Stimuli
stm_cnt = 1;
for iB = 1:2
    for iC = 1:numel(con)
        for iV = 1:numel(vow)
            
            if iB == 2 && (iC == 8 || iC == 12 || iC == 15)
                stm_cnt = stm_cnt + 1;
            else
                if iB == 1
                    vis_stm{stm_cnt} = [con{iC} vow{iV}];
                    if iV == 5; iswrd = cv_wrd_cut{find(strcmpi(cv_wrd_cut,[con{iC} 'A'])),2}; else iswrd = cv_wrd_cut{find(strcmpi(cv_wrd_cut,vis_stm{stm_cnt})),2}; end
                    if iswrd && iV == 5; wrd{end+1} = cv_wrd{find(strcmpi(cv_wrd_cut,[con{iC} 'A'])),2}; elseif iswrd; wrd{end+1} = cv_wrd_cut{find(strcmpi(cv_wrd_cut,vis_stm{stm_cnt})),2}; end
                elseif iB == 2
                    vis_stm{stm_cnt} = [vow{iV} con{iC}];
                    iswrd = vc_wrd_cut{find(strcmpi(vc_wrd_cut,vis_stm{stm_cnt})),2};
                    if iswrd; wrd{end+1} = vc_wrd{find(strcmpi(vc_wrd_cut,vis_stm{stm_cnt})),3}; end
                end
                
                con_con_num(stm_cnt) = iC;
                con_vow_num(stm_cnt) = iV;
                con_spk_num(stm_cnt) = spk_num_mtc(stm_cnt);
                con_bgr_ord(stm_cnt) = iB;
                con_wrd_num(stm_cnt) = iswrd;
                con_mes_con(stm_cnt) = mes_con(iC);
                con_mes_vow(stm_cnt) = mes_vow(iV);
                
                stm_cnt = stm_cnt + 1;
            end
        end
    end
end

con_con = ones(1,numel(vis_stm));

rmv_ind = find(cellfun(@isempty,vis_stm));

con_con(rmv_ind)     = []; 
vis_stm(rmv_ind)     = [];
con_con_num(rmv_ind) = [];
con_vow_num(rmv_ind) = [];
con_spk_num(rmv_ind) = [];
con_bgr_ord(rmv_ind) = [];
con_wrd_num(rmv_ind) = [];
con_mes_con(rmv_ind) = [];
con_mes_vow(rmv_ind) = [];

%% Create Auditory Stimuli
aud_stm = cell(1,numel(vis_stm));
for iW = 1:numel(vis_stm)
    
    if con_spk_num(iW) == 1
        if con_bgr_ord(iW) == 1 && con_vow_num(iW) == 5; aud_stm{iW} = ['F_' lower(vis_stm{iW}(1:end-1)) '.wav']; else aud_stm{iW} = ['F_' lower(vis_stm{iW})  '.wav']; end
%         aud_stm{iW} = ['F_' vis_stm{iW}  '.wav'];
    elseif con_spk_num(iW) == 2
        if con_bgr_ord(iW) == 1 && con_vow_num(iW) == 5; aud_stm{iW} = ['M_' lower(vis_stm{iW}(1:end-1)) '.wav']; else aud_stm{iW} = ['M_' lower(vis_stm{iW})  '.wav']; end
%         aud_stm{iW} = ['M_' vis_stm{iW}  '.wav'];
    end
    
end

%% Incongruous for SL-M

hld_vis_con_num = con_con_num;
hld_vis_vow_num = con_vow_num;
hld_aud_con_num = con_con_num;
hld_aud_vow_num = con_vow_num;

hld_spk_num = con_spk_num;

hld_vis_bgr_ord = con_bgr_ord;
hld_aud_bgr_ord = con_bgr_ord;

hld_vis_wrd_num = con_wrd_num;
hld_aud_wrd_num = con_wrd_num;

hld_vis_mes_con = con_mes_con;
hld_vis_mes_vow = con_mes_vow;
hld_aud_mes_con = con_mes_con;
hld_aud_mes_vow = con_mes_vow;

hld_vis_stm = vis_stm;
hld_aud_stm = aud_stm;
    
for iIN = 1:numel(aud_stm)
    
    inc_vis_stm{iIN} = hld_vis_stm{iIN};
    
    cmp = 0; fal = 1;
    while ~cmp
        
        inc_ind = floor((numel(hld_aud_stm)-1).*rand(1) + 1);
        
        if hld_vis_con_num(iIN) ~= hld_aud_con_num(inc_ind) && ...
                hld_vis_vow_num(iIN) ~= hld_aud_vow_num(inc_ind) && ...
                hld_vis_bgr_ord(iIN) == hld_aud_bgr_ord(inc_ind) && ...
                hld_vis_mes_con(iIN) ~= hld_aud_mes_con(inc_ind) && ...
                hld_vis_mes_vow(iIN) ~= hld_aud_mes_vow(inc_ind);
            cmp = 1;
        else
            fal = fal + 1; %             break;
            if fal > 50000
                fprintf(['FAIL COUNT:' num2str(fal) '\n'])
                break
            end
        end   
    end
                  
    fprintf('%i\n',iIN)
    
    con_inc(iIN) = 2; % congruent
    
    inc_vis_con_num(iIN) = hld_vis_con_num(iIN);
    inc_vis_vow_num(iIN) = hld_vis_vow_num(iIN);
    inc_aud_con_num(iIN) = hld_aud_con_num(inc_ind); hld_aud_con_num(inc_ind) = [];
    inc_aud_vow_num(iIN) = hld_aud_vow_num(inc_ind); hld_aud_vow_num(inc_ind) = [];
    
    inc_spk_num(iIN) = hld_spk_num(inc_ind); hld_spk_num(inc_ind) = [];
    
    inc_vis_bgr_ord(iIN) = con_bgr_ord(iIN);
    inc_aud_bgr_ord(iIN) = hld_aud_bgr_ord(inc_ind); hld_aud_bgr_ord(inc_ind) = [];
    
    inc_vis_wrd_num(iIN) = con_wrd_num(iIN);
    inc_aud_wrd_num(iIN) = hld_aud_wrd_num(inc_ind); hld_aud_wrd_num(inc_ind) = [];
    
    inc_vis_con_mes_con(iIN) = con_mes_con(iIN);
    inc_vis_con_mes_vow(iIN) = con_mes_vow(iIN);
    inc_aud_con_mes_con(iIN) = hld_aud_mes_con(inc_ind); hld_aud_mes_con(inc_ind) = [];
    inc_aud_con_mes_vow(iIN) = hld_aud_mes_vow(inc_ind); hld_aud_mes_vow(inc_ind) = [];
    
    inc_aud_stm{iIN} = hld_aud_stm{inc_ind}; hld_aud_stm(inc_ind) = [];

    
end

% sum(inc_vis_con_num==inc_aud_con_num)
% sum(inc_vis_vow_num==inc_aud_vow_num)
% sum(inc_vis_bgr_ord==inc_aud_bgr_ord)
% sum(inc_vis_con_mes_con==inc_aud_con_mes_con)
% sum(inc_vis_con_mes_vow==inc_aud_con_mes_vow)

% Put it all together

% Make Lists
SLM = [cat(1,vis_stm',inc_vis_stm'),cat(1,aud_stm',inc_aud_stm'),num2cell(cat(1,con_con',con_inc')),num2cell(cat(1,con_con_num',inc_vis_con_num')),num2cell(cat(1,con_vow_num',inc_vis_vow_num')), ...
    num2cell(cat(1,con_con_num',inc_aud_con_num')),num2cell(cat(1,con_vow_num',inc_aud_vow_num')),num2cell(cat(1,con_spk_num',inc_spk_num')),num2cell(cat(1,con_bgr_ord',inc_vis_bgr_ord')), ...
    num2cell(cat(1,con_wrd_num',inc_vis_wrd_num')),num2cell(cat(1,con_wrd_num',inc_aud_wrd_num')),num2cell(cat(1,con_mes_con',inc_vis_con_mes_con')), num2cell(cat(1,con_mes_vow',inc_vis_con_mes_vow')) ...
    num2cell(cat(1,con_mes_con',inc_aud_con_mes_con')),num2cell(cat(1,con_mes_vow',inc_aud_con_mes_vow'))];

% VisCon, VisVow, AudCon, AudVow, SpkNum, Bgr, VisWrd, AudWrd, VisConMesGrp, VisVowMesGrp, AudConMesGrp, VisVowMesGrp
stm_lst_rnd1 = ejk_createlist(SLM,[3 3 3 3 3 4 5 5 6 6 5 5 5],[3 4 5 6 7 8 9 10 11 12 13 14]);
ejk_presentationlist(stm_lst_rnd1,'/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_M/SL_M_Lists','SL_Auditory_First',4);

stm_lst_rnd2 = ejk_createlist(SLM,[3 3 3 3 3 4 5 5 6 6 5 5 5],[3 4 5 6 7 8 9 10 11 12 13 14]);
ejk_presentationlist(stm_lst_rnd2,'/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_M/SL_M_Lists','SL_Visual_First',4);

%% Copy Auditory Stimuli Over to Folder
auditory_stimuli = unique([stm_lst_rnd1(:,2);stm_lst_rnd2(:,2)]);

cv_auditory_stimuli = dir('/home/ekaestne/Downloads/CV_5_11_14/stimuli/*.wav');                         cv_auditory_stimuli = {cv_auditory_stimuli(:).name};
vc_auditory_stimuli = dir('/home/ekaestne/Downloads/VC_7_30_13/stimuli/AmpNorm_then_LengthNorm/*.wav'); vc_auditory_stimuli = {vc_auditory_stimuli(:).name};

for iAS = 1:numel(auditory_stimuli)
   cv_ind = find(strcmpi(auditory_stimuli{iAS},cv_auditory_stimuli));
   vc_ind = find(strcmpi(auditory_stimuli{iAS},vc_auditory_stimuli));
    
   if ~isempty(vc_ind)
       copyfile(['/home/ekaestne/Downloads/VC_7_30_13/stimuli/AmpNorm_then_LengthNorm/' auditory_stimuli{iAS}],'/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_Auditory_Stimuli')
   elseif ~isempty(cv_ind)
       copyfile(['/home/ekaestne/Downloads/CV_5_11_14/stimuli/' auditory_stimuli{iAS}],'/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_Auditory_Stimuli')
   end    
end

%% Incongruous for SL-L
% Choose the 3 groups
% Create VC
ttt11 = []; ttt21 = []; ttt31 = []; ttt32 = []; ttt12 = []; ttt22 = [];
ttt = repmat([1:9]',12,1);
ttt(:,2) = repmat([1:7 9:11 13:14]',9,1);
ttt(:,2) = sort(ttt(:,2));

ttt = reshape(ttt,9,12,2);

for iOT = 1:size(ttt,1)
    ttt(iOT,:,2) = circshift(squeeze(ttt(iOT,:,2))',iOT)';
    ttt11 = cat(1,ttt11,squeeze(ttt(iOT,1:3:size(ttt,2),:)));
    ttt21 = cat(1,ttt21,squeeze(ttt(iOT,2:3:size(ttt,2),:)));
    ttt31 = cat(1,ttt31,squeeze(ttt(iOT,3:3:size(ttt,2),:)));      
end
% 
% for iRM = 1:size(ttt11,1)
%     if (ttt11(iRM,2) == 8 || ttt11(iRM,2) == 12 || ttt11(iRM,2) == 15); rmv_ttt11(iRM) = 1; end
%     if (ttt21(iRM,2) == 8 || ttt21(iRM,2) == 12 || ttt21(iRM,2) == 15); rmv_ttt21(iRM) = 1; end
%     if (ttt31(iRM,2) == 8 || ttt31(iRM,2) == 12 || ttt31(iRM,2) == 15); rmv_ttt31(iRM) = 1; end
% end
% 
% ttt11(find(rmv_ttt11),:) = [];
% ttt21(find(rmv_ttt21),:) = [];
% ttt31(find(rmv_ttt31),:) = [];

% Put together measures
for iTT = 1:size(ttt11,1)
    fst_mis{iTT} = [vow{ttt11(iTT,1)} con{ttt11(iTT,2)}];
    fst_loc = find(strcmpi(fst_mis{iTT},vis_stm));
    fst_mis_con_num(iTT) = con_con_num(fst_loc);
    fst_mis_vow_num(iTT) = con_vow_num(fst_loc);
    fst_mis_spk_num(iTT) = con_spk_num(fst_loc);
    fst_mis_bgr_ord(iTT) = con_bgr_ord(fst_loc);
    fst_mis_wrd_num(iTT) = con_wrd_num(fst_loc);
    
    lst_mis{iTT} = [vow{ttt21(iTT,1)} con{ttt21(iTT,2)}];
    lst_loc = find(strcmpi(lst_mis{iTT},vis_stm));
    lst_mis_con_num(iTT) = con_con_num(lst_loc);
    lst_mis_vow_num(iTT) = con_vow_num(lst_loc);
    lst_mis_spk_num(iTT) = con_spk_num(lst_loc);
    lst_mis_bgr_ord(iTT) = con_bgr_ord(lst_loc);
    lst_mis_wrd_num(iTT) = con_wrd_num(lst_loc);
    
    bth_mis{iTT} = [vow{ttt31(iTT,1)} con{ttt31(iTT,2)}];
    bth_loc = find(strcmpi(bth_mis{iTT},vis_stm));
    bth_mis_con_num(iTT) = con_con_num(bth_loc);
    bth_mis_vow_num(iTT) = con_vow_num(bth_loc);
    bth_mis_spk_num(iTT) = con_spk_num(bth_loc);
    bth_mis_bgr_ord(iTT) = con_bgr_ord(bth_loc);
    bth_mis_wrd_num(iTT) = con_wrd_num(bth_loc);
end

% Create CV
ttt = repmat([1:15]',9,1);
ttt(:,2) = repmat([1:9]',15,1);
ttt(:,2) = sort(ttt(:,2));

ttt = reshape(ttt,15,9,2);

for iOT = 1:size(ttt,1)
    ttt(iOT,:,2)  = circshift(squeeze(ttt(iOT,:,2))',iOT)';
    ttt32         = cat(1,ttt32,squeeze(ttt(iOT,1:3:size(ttt,2),:)));
    ttt12         = cat(1,ttt12,squeeze(ttt(iOT,2:3:size(ttt,2),:)));
    ttt22         = cat(1,ttt22,squeeze(ttt(iOT,3:3:size(ttt,2),:)));      
end

% Put together measures
for iTT = 1:size(ttt12,1)
    fst_mis2{iTT} = [con{ttt12(iTT,1)} vow{ttt12(iTT,2)}];
    fst_loc2 = find(strcmpi(fst_mis2{iTT},vis_stm));
    fst_mis_con_num2(iTT) = con_con_num(fst_loc);
    fst_mis_vow_num2(iTT) = con_vow_num(fst_loc);
    fst_mis_spk_num2(iTT) = con_spk_num(fst_loc);
    fst_mis_bgr_ord2(iTT) = con_bgr_ord(fst_loc);
    fst_mis_wrd_num2(iTT) = con_wrd_num(fst_loc);
    
    lst_mis2{iTT} = [con{ttt22(iTT,1)} vow{ttt22(iTT,2)}];
    lst_loc2 = find(strcmpi(lst_mis2{iTT},vis_stm));
    lst_mis_con_num2(iTT) = con_con_num(lst_loc);
    lst_mis_vow_num2(iTT) = con_vow_num(lst_loc);
    lst_mis_spk_num2(iTT) = con_spk_num(lst_loc);
    lst_mis_bgr_ord2(iTT) = con_bgr_ord(lst_loc);
    lst_mis_wrd_num2(iTT) = con_wrd_num(lst_loc);
    
    bth_mis2{iTT} = [con{ttt32(iTT,1)} vow{ttt32(iTT,2)}];
    bth_loc2 = find(strcmpi(bth_mis2{iTT},vis_stm));
    bth_mis_con_num2(iTT) = con_con_num(bth_loc);
    bth_mis_vow_num2(iTT) = con_vow_num(bth_loc);
    bth_mis_spk_num2(iTT) = con_spk_num(bth_loc);
    bth_mis_bgr_ord2(iTT) = con_bgr_ord(bth_loc);
    bth_mis_wrd_num2(iTT) = con_wrd_num(bth_loc);
end

% Combine
fst_mis = cat(2,fst_mis,fst_mis2);
fst_mis_con_num = cat(2,fst_mis_con_num,fst_mis_con_num2);
fst_mis_vow_num = cat(2,fst_mis_vow_num,fst_mis_vow_num2);
fst_mis_spk_num = cat(2,fst_mis_spk_num,fst_mis_spk_num2);
fst_mis_bgr_ord = cat(2,fst_mis_bgr_ord,fst_mis_bgr_ord2);
fst_mis_wrd_num = cat(2,fst_mis_wrd_num,fst_mis_wrd_num2);

lst_mis = cat(2,lst_mis,lst_mis2);
lst_mis_con_num = cat(2,lst_mis_con_num,lst_mis_con_num2);
lst_mis_vow_num = cat(2,lst_mis_vow_num,lst_mis_vow_num2);
lst_mis_spk_num = cat(2,lst_mis_spk_num,lst_mis_spk_num2);
lst_mis_bgr_ord = cat(2,lst_mis_bgr_ord,lst_mis_bgr_ord2);
lst_mis_wrd_num = cat(2,lst_mis_wrd_num,lst_mis_wrd_num2);

bth_mis = cat(2,bth_mis,bth_mis2);
bth_mis_con_num = cat(2,bth_mis_con_num,bth_mis_con_num2);
bth_mis_vow_num = cat(2,bth_mis_vow_num,bth_mis_vow_num2);
bth_mis_spk_num = cat(2,bth_mis_spk_num,bth_mis_spk_num2);
bth_mis_bgr_ord = cat(2,bth_mis_bgr_ord,bth_mis_bgr_ord2);
bth_mis_wrd_num = cat(2,bth_mis_wrd_num,bth_mis_wrd_num2);

cell2csv('/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/VisStim/fst_mis',fst_mis')
cell2csv('/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/VisStim/lst_mis',lst_mis')
cell2csv('/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/VisStim/bth_mis',bth_mis')

% Mix up Sounds







% Split Auditory into Matching Stimuli
spt = ones(1,length(aud_stm)); spt(1:2:length(aud_stm)) = spt(1:2:length(aud_stm)) * 2;
spt = spt(randperm(numel(spt)));

aud_stm_hld = aud_stm;
aud_stm_fin = cell(1,numel(aud_stm));
vis_stm_fin = cell(1,numel(vis_stm));

rm_ind = zeros(1,numel(aud_stm));

for iIN = (numel(spt)/2)+1:numel(spt) % INcongruous 

    chs_whl = 1;
    while chs_whl
        pot_aud = find(~rm_ind); pot_aud = pot_aud(randperm(numel(pot_aud))); pot_aud = pot_aud(1);
        if ~strcmpi(aud_stm{pot_aud}(3:5),vis_stm(iIN));
           aud_stm_fin{iIN} = aud_stm{pot_aud}; rm_ind(aud_chs) = 1; chs_whl = 0;
        end
    end
    
    vis_stm_fin{iIN} = vis_stm{iIN};
        
    con_inc(iIN) = 2; % congruent
    vis_con(iIN) = con_num(iIN);% vis con
    vis_vow(iIN) = vow_num(iIN);% vis vow
    aud_con(iIN) = con_num(pot_aud);% aud con
    aud_vow(iIN) = vow_num(pot_aud);% aud vow
    aud_spk(iIN) = spk_num(pot_aud);% aud speaker
    vis_wrd(iIN) = wrd_num(iIN);% vis isword
    aud_wrd(iIN) = wrd_num(pot_aud);% aud isword
    vis_bgr(iIN) = bgr_ord(iIN); % vis bigram type
    aud_bgr(iIN) = bgr_ord(pot_aud); % aud bigram type
    
end    

% Save stimuli for visual inspection
stm_lst = [vis_stm_fin' aud_stm_fin' num2cell(con_inc)' num2cell(aud_spk)' num2cell(vis_wrd)' num2cell(aud_wrd)' num2cell(vis_bgr)' num2cell(vis_con)' num2cell(vis_vow)' num2cell(aud_bgr)' num2cell(aud_con)' num2cell(aud_vow)'];
cell2csv('/home/ekaestne/Desktop/Tasks/SL/NEWSTIM.csv',stm_lst)

% Test equal splits

% Make List
stm_lst_rnd = ejk_createlist(stm_lst,[3 3 5 5 4 2 2 4 2 2],[3 4 5 6 7 8 9 10 11 12]);
ejk_presentationlist(stm_lst_rnd,'/home/ekaestne/Desktop/Tasks/SL','NEWSTIM',3);









clear; clc;


%% Make lists for SL
con1 = {'N' 'T' 'S' 'W'}; % 
vow1 = {'EE' 'UH' 'OH' 'AY'}; %

con2 = {'M' 'D' 'H' 'F'}; %
vow2 = {'EE' 'OO' 'AH' 'UR'}; %

con3 = {'L' 'Y' 'R' 'K'}; %
vow3 = {'EH' 'AW' 'OY' 'IH'}; %

cv_wrd = mmil_readtext('/home/ekaestne/Desktop/Tasks/CV_words/CV_wordlist.csv');
cv_wrd_cut = cv_wrd(2:end,[1 3]);

%% Visual Stimuli
stm_cnt = 1;
for iC = 1:numel(con1)
    for iV = 1:numel(vow1)
        
        vis_stm1{stm_cnt} = [con1{iC} vow1{iV}];
        iswrd1 = cv_wrd_cut{find(strcmpi(cv_wrd_cut,vis_stm1{stm_cnt})),2};
        
        con_con_num1(stm_cnt) = iC;
        con_vow_num1(stm_cnt) = iV;
        con_wrd_num1(stm_cnt) = iswrd1;
        
        stm_cnt = stm_cnt + 1;
    end
end

stm_cnt = 1;
for iC = 1:numel(con2)
    for iV = 1:numel(vow2)
        
        vis_stm2{stm_cnt} = [con2{iC} vow2{iV}];
        iswrd2 = cv_wrd_cut{find(strcmpi(cv_wrd_cut,vis_stm2{stm_cnt})),2};
        
        con_con_num2(stm_cnt) = iC+4;
        con_vow_num2(stm_cnt) = iV+4;
        con_wrd_num2(stm_cnt) = iswrd2;
        
        stm_cnt = stm_cnt + 1;
    end
end

stm_cnt = 1;
for iC = 1:numel(con3)
    for iV = 1:numel(vow3)
        
        vis_stm3{stm_cnt} = [con3{iC} vow3{iV}];
        iswrd3 = cv_wrd_cut{find(strcmpi(cv_wrd_cut,vis_stm3{stm_cnt})),2};
        
        con_con_num3(stm_cnt) = iC+8;
        con_vow_num3(stm_cnt) = iV+8;
        con_wrd_num3(stm_cnt) = iswrd3;
        
        stm_cnt = stm_cnt + 1;
    end
end

%% Create Auditory Stimuli
aud_stm = cell(1,numel(vis_stm1));
for iW = 1:numel(vis_stm1)
    
    aud_stm1{iW} = ['F_' lower(vis_stm1{iW})  '.wav'];
    
    if isempty(strfind(vis_stm2{iW},'AH'))
        aud_stm2{iW} = ['F_' lower(vis_stm2{iW})  '.wav'];
    else
        aud_stm2{iW} = ['F_' lower(vis_stm2{iW}(1:end-1))  '.wav'];
    end
    
    aud_stm_nse1{iW} = ['F_' lower(vis_stm1{iW})  '_matched_noise_TB20_Fc50-5000.wav'];
    
    if isempty(strfind(vis_stm2{iW},'AH'))
    aud_stm_nse2{iW} = ['F_' lower(vis_stm2{iW})  '_matched_noise_TB20_Fc50-5000.wav'];
    else
        aud_stm_nse2{iW} = ['F_' lower(vis_stm2{iW}(1:end-1))  '_matched_noise_TB20_Fc50-5000.wav'];
    end
end
for iW = 1:numel(vis_stm3)
    aud_stm3{iW} = ['F_' lower(vis_stm3{iW})  '.wav'];
    aud_stm_nse3{iW} = ['F_' lower(vis_stm3{iW})  '_matched_noise_TB20_Fc50-5000.wav'];
end
%% Copy Auditory Stimuli Over to Folder
% cv_auditory_stimuli = dir('/home/ekaestne/Downloads/CV_5_11_14/stimuli/*.wav'); cv_auditory_stimuli = {cv_auditory_stimuli(:).name};
% 
% for iAS = 1:numel(aud_stm1)
%     cv_ind1 = find(strcmpi(aud_stm1{iAS},cv_auditory_stimuli));
%     cv_ind2 = find(strcmpi(aud_stm2{iAS},cv_auditory_stimuli));
%     
%     cv_ind_nse1 = find(strcmpi(aud_stm_nse1{iAS},cv_auditory_stimuli));
%     cv_ind_nse2 = find(strcmpi(aud_stm_nse2{iAS},cv_auditory_stimuli));
%         
%     copyfile(['/home/ekaestne/Downloads/CV_5_11_14/stimuli/' cv_auditory_stimuli{cv_ind1}],'/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_F_Stimuli')
%     copyfile(['/home/ekaestne/Downloads/CV_5_11_14/stimuli/' cv_auditory_stimuli{cv_ind2}],'/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_F_Stimuli')
%         
%     copyfile(['/home/ekaestne/Downloads/CV_5_11_14/stimuli/' cv_auditory_stimuli{cv_ind_nse1}],'/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_F_Stimuli')
%     copyfile(['/home/ekaestne/Downloads/CV_5_11_14/stimuli/' cv_auditory_stimuli{cv_ind_nse2}],'/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_F_Stimuli')
% end
% 
% for iAS = 1:numel(aud_stm3)
%     cv_ind3 = find(strcmpi(aud_stm3{iAS},cv_auditory_stimuli));
%     cv_ind_nse3 = find(strcmpi(aud_stm_nse3{iAS},cv_auditory_stimuli));
%     copyfile(['/home/ekaestne/Downloads/CV_5_11_14/stimuli/' cv_auditory_stimuli{cv_ind3}],'/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_F_Stimuli')
%     copyfile(['/home/ekaestne/Downloads/CV_5_11_14/stimuli/' cv_auditory_stimuli{cv_ind_nse3}],'/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_F_Stimuli')
% end

 %% Incongruous for SL-F, list 1
out_cmp = 0;

for iN = 1:4
    while ~out_cmp
        
        hld_vis_con_num = con_con_num1;
        hld_vis_vow_num = con_vow_num1;
        hld_aud_con_num = con_con_num1;
        hld_aud_vow_num = con_vow_num1;
        
        hld_vis_wrd_num = con_wrd_num1;
        hld_aud_wrd_num = con_wrd_num1;
        
        hld_vis_stm = vis_stm1;
        hld_aud_stm = aud_stm1;
        
        for iIN = 1:numel(vis_stm1)
            
            inc_vis_stm{iIN} = hld_vis_stm{iIN};
            
            inn_cmp = 0; fal = 1;
            while ~inn_cmp
                
                inc_ind = floor((numel(hld_aud_stm)).*rand(1) + 1);
                
                if hld_vis_con_num(iIN) ~= hld_aud_con_num(inc_ind) && ...
                        hld_vis_vow_num(iIN) ~= hld_aud_vow_num(inc_ind)
                    inn_cmp = 1;
                else
                    fal = fal + 1; %             break;
                    if fal > 50000
                        fprintf(['FAIL COUNT:' num2str(fal) '\n'])
                        break
                    end
                end
            end
            
            %     fprintf('%i\n',iIN)
            
            if fal < 50000;
                con_inc(iIN) = 2; % congruent
                
                inc_vis_con_num(iIN) = hld_vis_con_num(iIN);
                inc_vis_vow_num(iIN) = hld_vis_vow_num(iIN);
                inc_aud_con_num(iIN) = hld_aud_con_num(inc_ind); hld_aud_con_num(inc_ind) = [];
                inc_aud_vow_num(iIN) = hld_aud_vow_num(inc_ind); hld_aud_vow_num(inc_ind) = [];
                
                inc_vis_wrd_num(iIN) = hld_vis_wrd_num(iIN);
                inc_aud_wrd_num(iIN) = hld_aud_wrd_num(inc_ind); hld_aud_wrd_num(inc_ind) = [];
                
                inc_aud_stm{iIN} = hld_aud_stm{inc_ind}; hld_aud_stm(inc_ind) = [];
                if isempty(hld_aud_con_num); out_cmp = 1; end
            else
                break
            end
            
        end
    end
    
    lst1.(['num_' num2str(iN)]) = [vis_stm1' aud_stm1' num2cell(ones(numel(vis_stm1),1)*1) num2cell(inc_vis_con_num') num2cell(inc_vis_vow_num') num2cell(con_con_num1') num2cell(con_vow_num1'); ...
        vis_stm1' inc_aud_stm' num2cell(ones(numel(vis_stm1),1)*2) num2cell(inc_vis_con_num') num2cell(inc_vis_vow_num') num2cell(inc_aud_con_num') num2cell(inc_aud_vow_num'); ...
        vis_stm1' aud_stm1' num2cell(ones(numel(vis_stm1),1)*3) num2cell(inc_vis_con_num') num2cell(inc_vis_vow_num') num2cell(con_con_num1') num2cell(con_vow_num1'); ...
        vis_stm1' aud_stm_nse1' num2cell(ones(numel(vis_stm1),1)*4) num2cell(inc_vis_con_num') num2cell(inc_vis_vow_num') num2cell(con_con_num1') num2cell(con_vow_num1')];
    
end

%% Incongruous for SL-F, list 2
out_cmp = 0;

for iN = 1:4
    while ~out_cmp
        
        hld_vis_con_num = con_con_num2;
        hld_vis_vow_num = con_vow_num2;
        hld_aud_con_num = con_con_num2;
        hld_aud_vow_num = con_vow_num2;
        
        hld_vis_wrd_num = con_wrd_num2;
        hld_aud_wrd_num = con_wrd_num2;
        
        hld_vis_stm = vis_stm2;
        hld_aud_stm = aud_stm2;
        
        for iIN = 1:numel(vis_stm2)
            
            inc_vis_stm{iIN} = hld_vis_stm{iIN};
            
            inn_cmp = 0; fal = 1;
            while ~inn_cmp
                
                inc_ind = floor((numel(hld_aud_stm)).*rand(1) + 1);
                
                if hld_vis_con_num(iIN) ~= hld_aud_con_num(inc_ind) && ...
                        hld_vis_vow_num(iIN) ~= hld_aud_vow_num(inc_ind)
                    inn_cmp = 1;
                else
                    fal = fal + 1; %             break;
                    if fal > 50000
                        fprintf(['FAIL COUNT:' num2str(fal) '\n'])
                        break
                    end
                end
            end
            
            %     fprintf('%i\n',iIN)
            
            if fal < 50000;
                con_inc(iIN) = 2; % congruent
                
                inc_vis_con_num(iIN) = hld_vis_con_num(iIN);
                inc_vis_vow_num(iIN) = hld_vis_vow_num(iIN);
                inc_aud_con_num(iIN) = hld_aud_con_num(inc_ind); hld_aud_con_num(inc_ind) = [];
                inc_aud_vow_num(iIN) = hld_aud_vow_num(inc_ind); hld_aud_vow_num(inc_ind) = [];
                
                inc_vis_wrd_num(iIN) = hld_vis_wrd_num(iIN);
                inc_aud_wrd_num(iIN) = hld_aud_wrd_num(inc_ind); hld_aud_wrd_num(inc_ind) = [];
                
                inc_aud_stm{iIN} = hld_aud_stm{inc_ind}; hld_aud_stm(inc_ind) = [];
                if isempty(hld_aud_con_num); out_cmp = 1; end
            else
                break
            end
            
        end
    end
    
    lst2.(['num_' num2str(iN)]) = [vis_stm2' aud_stm2' num2cell(ones(numel(vis_stm2),1)*1) num2cell(inc_vis_con_num') num2cell(inc_vis_vow_num') num2cell(con_con_num2') num2cell(con_vow_num2'); ...
        vis_stm2' inc_aud_stm' num2cell(ones(numel(vis_stm2),1)*2) num2cell(inc_vis_con_num') num2cell(inc_vis_vow_num') num2cell(inc_aud_con_num') num2cell(inc_aud_vow_num'); ...
        vis_stm2' aud_stm2' num2cell(ones(numel(vis_stm2),1)*3) num2cell(inc_vis_con_num') num2cell(inc_vis_vow_num') num2cell(con_con_num2') num2cell(con_vow_num2'); ...
        vis_stm2' aud_stm_nse2' num2cell(ones(numel(vis_stm2),1)*4) num2cell(inc_vis_con_num') num2cell(inc_vis_vow_num') num2cell(con_con_num2') num2cell(con_vow_num2')];
    
end

%% Incongruous for SL-F, list 3
out_cmp = 0;

clear inc_vis_con_num inc_vis_vow_num inc_aud_con_num inc_aud_vow_num inc_vis_wrd_num inc_aud_wrd_num inc_aud_stm

for iN = 1:4
    while ~out_cmp
        
        hld_vis_con_num = con_con_num3;
        hld_vis_vow_num = con_vow_num3;
        hld_aud_con_num = con_con_num3;
        hld_aud_vow_num = con_vow_num3;
        
        hld_vis_wrd_num = con_wrd_num3;
        hld_aud_wrd_num = con_wrd_num3;
        
        hld_vis_stm = vis_stm3;
        hld_aud_stm = aud_stm3;
        
        for iIN = 1:numel(vis_stm3)
            
            inc_vis_stm{iIN} = hld_vis_stm{iIN};
            
            inn_cmp = 0; fal = 1;
            while ~inn_cmp
                
                inc_ind = floor((numel(hld_aud_stm)).*rand(1) + 1);
                
                if hld_vis_con_num(iIN) ~= hld_aud_con_num(inc_ind) && ...
                        hld_vis_vow_num(iIN) ~= hld_aud_vow_num(inc_ind)
                    inn_cmp = 1;
                else
                    fal = fal + 1; %             break;
                    if fal > 50000
                        fprintf(['FAIL COUNT:' num2str(fal) '\n'])
                        break
                    end
                end
            end
            
            %     fprintf('%i\n',iIN)
            
            if fal < 50000;
                con_inc(iIN) = 2; % congruent
                
                inc_vis_con_num(iIN) = hld_vis_con_num(iIN);
                inc_vis_vow_num(iIN) = hld_vis_vow_num(iIN);
                inc_aud_con_num(iIN) = hld_aud_con_num(inc_ind); hld_aud_con_num(inc_ind) = [];
                inc_aud_vow_num(iIN) = hld_aud_vow_num(inc_ind); hld_aud_vow_num(inc_ind) = [];
                
                inc_vis_wrd_num(iIN) = hld_vis_wrd_num(iIN);
                inc_aud_wrd_num(iIN) = hld_aud_wrd_num(inc_ind); hld_aud_wrd_num(inc_ind) = [];
                
                inc_aud_stm{iIN} = hld_aud_stm{inc_ind}; hld_aud_stm(inc_ind) = [];
                if isempty(hld_aud_con_num); out_cmp = 1; end
            else
                break
            end
            
        end
    end
    
    lst3.(['num_' num2str(iN)]) = [vis_stm3' aud_stm3' num2cell(ones(numel(vis_stm3),1)*1) num2cell(inc_vis_con_num') num2cell(inc_vis_vow_num') num2cell(con_con_num3') num2cell(con_vow_num3'); ...
        vis_stm3' inc_aud_stm' num2cell(ones(numel(vis_stm3),1)*2) num2cell(inc_vis_con_num') num2cell(inc_vis_vow_num') num2cell(inc_aud_con_num') num2cell(inc_aud_vow_num'); ...
        vis_stm3' aud_stm3' num2cell(ones(numel(vis_stm3),1)*3) num2cell(inc_vis_con_num') num2cell(inc_vis_vow_num') num2cell(con_con_num3') num2cell(con_vow_num3'); ...
        vis_stm3' aud_stm_nse3' num2cell(ones(numel(vis_stm3),1)*4) num2cell(inc_vis_con_num') num2cell(inc_vis_vow_num') num2cell(con_con_num3') num2cell(con_vow_num3')];
    
end

%% Put it all together
lst1_1 = [lst1.num_1 ; lst1.num_2];
lst1_1 = ejk_createlist(lst1_1,[2 2 2 2 2],[3:size(lst1_1,2)]);
ejk_presentationlist(lst1_1,'/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_F_Lists','SL_Vis1st_1_1',1);

lst1_2 = [lst1.num_3 ; lst1.num_4];
lst1_2 = ejk_createlist(lst1_2,[2 2 2 2 2],[3:size(lst1_2,2)]);
ejk_presentationlist(lst1_2,'/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_F_Lists','SL_Vis1st_1_2',1);

lst2_1 = [lst2.num_1 ; lst2.num_2];
lst2_1 = ejk_createlist(lst2_1,[2 2 2 2 2],[3:size(lst2_1,2)]);
ejk_presentationlist(lst2_1,'/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_F_Lists','SL_Vis1st_2_1',1);

lst2_2 = [lst2.num_3 ; lst2.num_4];
lst2_2 = ejk_createlist(lst2_2,[2 2 2 2 2],[3:size(lst2_2,2)]);
ejk_presentationlist(lst2_2,'/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_F_Lists','SL_Vis1st_2_2',1);

lst3_1 = [lst3.num_1 ; lst3.num_2];
lst3_1 = ejk_createlist(lst3_1,[2 2 2 2 2],[3:size(lst3_1,2)]);
ejk_presentationlist(lst3_1,'/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_F_Lists','SL_Vis1st_2_1',1);

lst3_2 = [lst2.num_3 ; lst2.num_4];
lst3_2 = ejk_createlist(lst3_2,[2 2 2 2 2],[3:size(lst3_2,2)]);
ejk_presentationlist(lst3_2,'/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_F_Lists','SL_Vis1st_2_2',1);

%% Make Practice Run
clear; clc;

tst_con = {'M' 'B' 'R'};
tst_vow = {'EE' 'UR' 'OH'};

cv_wrd = mmil_readtext('/home/ekaestne/Desktop/Tasks/CV_words/CV_wordlist.csv');
cv_wrd_cut = cv_wrd(2:end,[1 3]);

stm_cnt = 1;
for iC = 1:numel(tst_con)
    for iV = 1:numel(tst_vow)
        
        tst_stm{stm_cnt} = [tst_con{iC} tst_vow{iV}];
        iswrd1 = cv_wrd_cut{find(strcmpi(cv_wrd_cut,tst_stm{stm_cnt})),2};
        
        tst_con_num1(stm_cnt) = iC;
        tst_vow_num1(stm_cnt) = iV;
        tst_wrd_num1(stm_cnt) = iswrd1;
        
        stm_cnt = stm_cnt + 1;
    end
end

aud_stm = cell(1,numel(tst_stm));
for iW = 1:numel(tst_stm)
    
    tst_aud_stm1{iW} = ['F_' lower(tst_stm{iW})  '.wav'];
    
    tst_stm_nse1{iW} = ['F_' lower(tst_stm{iW})  '_matched_noise_TB20_Fc50-5000.wav'];
    
end

out_cmp = 0;

while ~out_cmp
    
    hld_vis_con_num = tst_con_num1;
    hld_vis_vow_num = tst_vow_num1;
    hld_aud_con_num = tst_con_num1;
    hld_aud_vow_num = tst_vow_num1;
    
    hld_vis_wrd_num = tst_wrd_num1;
    hld_aud_wrd_num = tst_wrd_num1;
    
    hld_vis_stm = tst_stm;
    hld_aud_stm = tst_aud_stm1;
    
    for iIN = 1:numel(tst_stm)
        
        inc_vis_stm{iIN} = hld_vis_stm{iIN};
        
        inn_cmp = 0; fal = 1;
        while ~inn_cmp
            
            inc_ind = floor((numel(hld_aud_stm)).*rand(1) + 1);
            
            if hld_vis_con_num(iIN) ~= hld_aud_con_num(inc_ind) && ...
                    hld_vis_vow_num(iIN) ~= hld_aud_vow_num(inc_ind)
                inn_cmp = 1;
            else
                fal = fal + 1; %             break;
                if fal > 50000
                    fprintf(['FAIL COUNT:' num2str(fal) '\n'])
                    break
                end
            end
        end
        
        %     fprintf('%i\n',iIN)
        
        if fal < 50000;
            con_inc(iIN) = 2; % congruent
            
            inc_vis_con_num(iIN) = hld_vis_con_num(iIN);
            inc_vis_vow_num(iIN) = hld_vis_vow_num(iIN);
            inc_aud_con_num(iIN) = hld_aud_con_num(inc_ind); hld_aud_con_num(inc_ind) = [];
            inc_aud_vow_num(iIN) = hld_aud_vow_num(inc_ind); hld_aud_vow_num(inc_ind) = [];
            
            inc_vis_wrd_num(iIN) = hld_vis_wrd_num(iIN);
            inc_aud_wrd_num(iIN) = hld_aud_wrd_num(inc_ind); hld_aud_wrd_num(inc_ind) = [];
            
            inc_aud_stm{iIN} = hld_aud_stm{inc_ind}; hld_aud_stm(inc_ind) = [];
            if isempty(hld_aud_con_num); out_cmp = 1; end
        else
            break
        end
        
    end
end

tst_lst = [tst_stm' tst_aud_stm1' num2cell(ones(numel(tst_stm),1)*1) num2cell(inc_vis_con_num') num2cell(inc_vis_vow_num') num2cell(tst_con_num1') num2cell(tst_vow_num1'); ...
    tst_stm' inc_aud_stm' num2cell(ones(numel(tst_stm),1)*2) num2cell(inc_vis_con_num') num2cell(inc_vis_vow_num') num2cell(inc_aud_con_num') num2cell(inc_aud_vow_num'); ...
    tst_stm' tst_aud_stm1' num2cell(ones(numel(tst_stm),1)*3) num2cell(inc_vis_con_num') num2cell(inc_vis_vow_num') num2cell(tst_con_num1') num2cell(tst_vow_num1'); ...
    tst_stm' tst_stm_nse1' num2cell(ones(numel(tst_stm),1)*4) num2cell(inc_vis_con_num') num2cell(inc_vis_vow_num') num2cell(tst_con_num1') num2cell(tst_vow_num1')];

tst_lst = ejk_createlist(tst_lst,[2 2 2 2 2],[3:size(tst_lst,2)]);
ejk_presentationlist(tst_lst,'/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_F_Lists','tst_lst',1);

%% Make Pronounciation Run
clear; clc;

% T1-B2 OW2- M1-N2

con1 = {'N' 'B' 'S' 'W'}; % 
vow1 = {'IH' 'UH' 'OH' 'AY'}; %

con2 = {'M' 'T' 'H' 'F'}; %
vow2 = {'EE' 'OO' 'AH' 'UR'}; %

con3 = {'L' 'Y' 'R' 'K'}; %
vow3 = {'EH' 'AW' 'OY'}; %

vis_lst = [strcat(con2,vow1) ...
    strcat(con3,vow1) ...
    strcat(con1,vow2) ...
    strcat(con3,vow2) ...
    strcat(con1(1:3),vow3) ...
    strcat(con2(1:3),vow3)];
aud_lst = strcat('F_',vis_lst,'.wav');

cell2csv('/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_F_Lists/Prounounciation_Practice.csv',[vis_lst',aud_lst'])

%% WordQ
cv_wrd = mmil_readtext('/home/ekaestne/Desktop/Tasks/CV_words/CV_wordlist.csv');
cv_wrd_cut = cv_wrd(2:end,[1 3]);

stm_cnt = 1;
for iC = 1:numel(con1)
    for iV = 1:numel(vow1)
        
        vis_stm1{stm_cnt} = [con1{iC} vow1{iV}];
        iswrd1 = cv_wrd_cut{find(strcmpi(cv_wrd_cut,vis_stm1{stm_cnt})),2};
        vis_stm2{stm_cnt} = [con2{iC} vow1{iV}];
        iswrd2 = cv_wrd_cut{find(strcmpi(cv_wrd_cut,vis_stm2{stm_cnt})),2};
        vis_stm3{stm_cnt} = [con3{iC} vow1{iV}];
        iswrd3 = cv_wrd_cut{find(strcmpi(cv_wrd_cut,vis_stm3{stm_cnt})),2};
        
        con_wrd_num1(stm_cnt) = iswrd1;
        con_wrd_num2(stm_cnt) = iswrd2;
        con_wrd_num3(stm_cnt) = iswrd3;
        
        stm_cnt = stm_cnt + 1;
    end
end
[vis_stm1' num2cell(con_wrd_num1)']% vis_stm2' num2cell(con_wrd_num2)' vis_stm3' num2cell(con_wrd_num3)']


stm_cnt = 1;
for iC = 1:numel(con2)
    for iV = 1:numel(vow2)
        
        vis_stm1{stm_cnt} = [con1{iC} vow2{iV}];
        iswrd1 = cv_wrd_cut{find(strcmpi(cv_wrd_cut,vis_stm1{stm_cnt})),2};
        vis_stm2{stm_cnt} = [con2{iC} vow2{iV}];
        iswrd2 = cv_wrd_cut{find(strcmpi(cv_wrd_cut,vis_stm2{stm_cnt})),2};
        vis_stm3{stm_cnt} = [con3{iC} vow2{iV}];
        iswrd3 = cv_wrd_cut{find(strcmpi(cv_wrd_cut,vis_stm3{stm_cnt})),2};
        
        con_wrd_num1(stm_cnt) = iswrd1;
        con_wrd_num2(stm_cnt) = iswrd2;
        con_wrd_num3(stm_cnt) = iswrd3;
        
        stm_cnt = stm_cnt + 1;
    end
end
[vis_stm2' num2cell(con_wrd_num2)'] %[vis_stm1' num2cell(con_wrd_num1)' vis_stm2' num2cell(con_wrd_num2)' vis_stm3' num2cell(con_wrd_num3)']
clear vis_stm1 vis_stm2 vis_stm3 con_wrd_num1 con_wrd_num2 con_wrd_num3

stm_cnt = 1;
for iC = 1:numel(con3)
    for iV = 1:numel(vow3)
        
        vis_stm1{stm_cnt} = [con1{iC} vow3{iV}];
        iswrd1 = cv_wrd_cut{find(strcmpi(cv_wrd_cut,vis_stm1{stm_cnt})),2};
        vis_stm2{stm_cnt} = [con2{iC} vow3{iV}];
        iswrd2 = cv_wrd_cut{find(strcmpi(cv_wrd_cut,vis_stm2{stm_cnt})),2};
        vis_stm3{stm_cnt} = [con3{iC} vow3{iV}];
        iswrd3 = cv_wrd_cut{find(strcmpi(cv_wrd_cut,vis_stm3{stm_cnt})),2};
        
        con_wrd_num1(stm_cnt) = iswrd1;
        con_wrd_num2(stm_cnt) = iswrd2;
        con_wrd_num3(stm_cnt) = iswrd3;
        
        stm_cnt = stm_cnt + 1;
    end
end
[vis_stm3' num2cell(con_wrd_num3)'] %[vis_stm1' num2cell(con_wrd_num1)' vis_stm2' num2cell(con_wrd_num2)' vis_stm3' num2cell(con_wrd_num3)']

tot_vow = [vow1 vow2 vow3];
tot_con = [con1 con2 con3];

stm_cnt = 1;
for iC = 1:numel(tot_vow)
    for iV = 1:numel(tot_con)
        
        vis_stm1{stm_cnt} = [tot_con{iC} tot_vow{iV}];
        iswrd1 = cv_wrd_cut{find(strcmpi(cv_wrd_cut,vis_stm1{stm_cnt})),2};
        
        con_wrd_num1(stm_cnt) = iswrd1;
        
        stm_cnt = stm_cnt + 1;
    end
end

[vis_stm1' num2cell(con_wrd_num1)']






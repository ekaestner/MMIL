con1 = {'N'  'T'  'S'  'W'}; % 
vow1 = {'EE' 'UH' 'OH' 'AY'}; %

cnt = 1;
for iC = 1:numel(con1)
    for iV = 1:numel(vow1)
        stm1{cnt} = [con1{iC} vow1{iV}];
        cnt = cnt+1;
    end
end

con2 = {'M'  'D'  'H'  'F'}; %
vow2 = {'EE' 'OO' 'AH' 'UR'}; %

cnt = 1;
for iC = 1:numel(con2)
    for iV = 1:numel(vow2)
        stm2{cnt} = [con2{iC} vow2{iV}];
        cnt = cnt+1;
    end
end

con3 = {'L'  'Y'  'R'  'K'}; %
vow3 = {'EH' 'AW' 'OY' 'IH'}; %

cnt = 1;
for iC = 1:numel(con3)
    for iV = 1:numel(vow3)
        stm3{cnt} = [con3{iC} vow3{iV}];
        cnt = cnt+1;
    end
end

stm = [stm1 stm2 stm3];

%% Copy Over Stimuli
for iCP = 1:numel(stm)
    copyfile(['/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_F_Stimuliv2/' 'F_' lower(stm{iCP}) '.wav'],['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Stimuli/Female/Voices' '/' 'F_' lower(stm{iCP}) '.wav'])
    copyfile(['/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_F_Stimuliv2/noise/' 'F_' lower(stm{iCP}) '_matched_noise_TB12_Fc20-20000_nobroadband_log.wav'],['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Stimuli/Female/Noise' '/' 'F_' lower(stm{iCP}) '_matched_noise_TB12_Fc20-20000_nobroadband_log.wav'])
    
%     copyfile(['/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_F_Stimuliv2' 'M_' lower(stm{iCP}) '.wav'],['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim/Stimuli/Male/Voices' '/' 'M_' lower(stm{iCP}) '.wav'])
%     copyfile(['/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_F_Stimuliv2' 'M_' lower(stm{iCP}) '_matched_noise_TB20_Fc50-5000.wav'],['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim/Stimuli/Male/Noise' '/' 'M_' lower(stm{iCP}) '_matched_noise_TB20_Fc50-5000.wav'])
end

%% Copy Over Figures
% for iCP = 1:numel(stm)
%    copyfile(['/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/BG/Stimuli/Figures' '/' 'F_' lower(stm{iCP}) '.png'],['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Figures/Female/Voices' '/' 'F_' lower(stm{iCP}) '.png']);
% %    copyfile(['/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/BG/Stimuli/Figures' '/' 'M_' lower(stm{iCP}) '.png'],['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Figures/Male/Voices' '/' 'M_' lower(stm{iCP}) '.png']);
% end

%% Copy Over Timings
for iCP = 1:numel(stm)
   copyfile(['/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/BG/Stimuli/Output' '/' 'F_' lower(stm{iCP}) '.txt'],['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Timings/Female/Voices' '/' 'F_' lower(stm{iCP}) '.txt']);
   copyfile(['/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/BG/Stimuli/Output' '/' 'M_' lower(stm{iCP}) '.txt'],['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Timings/Male/Voices' '/' 'M_' lower(stm{iCP}) '.txt']);
end

%% Make Timings
cmu_dic = mmil_readtext(['/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/BG/Stimuli' '/' 'CMU_Dictionary_Addon'],[' ']); cmu_dic(:,2) = [];
f_tme_hld = cell(numel(stm),3);
m_tme_hld = cell(numel(stm),3);
for iOF = 1:numel(stm)
   
    f_txt_tmp = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Timings/Female/Voices' '/' 'F_' lower(stm{iOF}) '.txt']);
    m_txt_tmp = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Timings/Male/Voices' '/' 'M_' lower(stm{iOF}) '.txt']);
   
    if strcmpi(f_txt_tmp{end,:},'sp'); f_dic_ind = find(strcmpi(cmu_dic(:,1),f_txt_tmp{end-3,:})); else f_dic_ind = find(strcmpi(cmu_dic(:,1),f_txt_tmp{end,:})); end  
    if strcmpi(m_txt_tmp{end,:},'sp'); m_dic_ind = find(strcmpi(cmu_dic(:,1),m_txt_tmp{end-3,:})); else m_dic_ind = find(strcmpi(cmu_dic(:,1),m_txt_tmp{end,:})); end
    
    if f_dic_ind ~= m_dic_ind; error('f & m do not match'); end
    
    f_fst_phn = cmu_dic(f_dic_ind,2);
    f_scd_phn = cmu_dic(f_dic_ind,3);
        
    m_fst_phn = cmu_dic(m_dic_ind,2);
    m_scd_phn = cmu_dic(m_dic_ind,3);
    
    f_fst_phn_lne = find(strcmpi(f_txt_tmp,f_fst_phn));
    f_scd_phn_lne = find(strcmpi(f_txt_tmp,f_scd_phn));
    
    if isempty(f_fst_phn_lne)
        f_fst_phn_lne = f_scd_phn_lne - 3;
    elseif isempty(f_scd_phn_lne)
        f_scd_phn_lne = f_fst_phn_lne + 3;
    end
    
    m_fst_phn_lne = find(strcmpi(m_txt_tmp,m_fst_phn));
    m_scd_phn_lne = find(strcmpi(m_txt_tmp,m_scd_phn));
    
    if isempty(m_fst_phn_lne)
        m_fst_phn_lne = m_scd_phn_lne - 3;
    elseif isempty(m_scd_phn_lne)
        m_scd_phn_lne = m_fst_phn_lne + 3;
    end
    
    if ~isempty(f_fst_phn_lne) && ~isempty(f_scd_phn_lne)
        f_fst_phn_tme = f_txt_tmp(f_fst_phn_lne-2);
        f_scd_phn_tme = f_txt_tmp(f_scd_phn_lne-2);
        f_tme_hld(iOF,:) = [stm{iOF} f_fst_phn_tme f_scd_phn_tme];
    else
        
    end
    
    if ~isempty(m_fst_phn_lne) && ~isempty(m_scd_phn_lne)
        m_fst_phn_tme = m_txt_tmp(m_fst_phn_lne-2);
        m_scd_phn_tme = m_txt_tmp(m_scd_phn_lne-2);
        m_tme_hld(iOF,:) = [stm{iOF} m_fst_phn_tme m_scd_phn_tme];
    else
        
    end
    
end

cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/female_phoneme_timing.csv',f_tme_hld)
cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/male_phoneme_timing.csv',m_tme_hld)

%% Make Plots
f_tme_hld = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/female_phoneme_timing.csv'); [~,ind] = sort(f_tme_hld(:,1)); f_tme_hld = f_tme_hld(ind,:);
% m_tme_hld = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/male_phoneme_timing.csv');   [~,ind] = sort(m_tme_hld(:,1)); m_tme_hld = m_tme_hld(ind,:);

f_vce_wav_fle = dir(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Stimuli' '/' 'Female' '/' 'Voices' '/' '*.wav']); f_vce_wav_fle = {f_vce_wav_fle(:).name};
f_nse_wav_fle = dir(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Stimuli' '/' 'Female' '/' 'Noise' '/' '*.wav']);  f_nse_wav_fle = {f_nse_wav_fle(:).name};
% m_vce_wav_fle = dir(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Stimuli' '/' 'Male' '/' 'Voices' '/' '*.wav']);   m_vce_wav_fle = {m_vce_wav_fle(:).name};
% m_nse_wav_fle = dir(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Stimuli' '/' 'Male' '/' 'Noise' '/' '*.wav']);    m_nse_wav_fle = {m_nse_wav_fle(:).name};

for iFG = 1:numel(f_vce_wav_fle)
    
    [f_vce_wav,fs_f_vce_wav] = audioread(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Stimuli' '/' 'Female' '/' 'Voices' '/' f_vce_wav_fle{iFG}]);
    [f_nse_wav,fs_f_nse_wav] = audioread(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Stimuli' '/' 'Female' '/' 'Noise' '/' f_nse_wav_fle{iFG}]);
%     [m_vce_wav,fs_m_vce_wav] = audioread(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Stimuli' '/' 'Male' '/' 'Voices' '/' m_vce_wav_fle{iFG}]);
%     [m_nse_wav,fs_m_nse_wav] = audioread(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Stimuli' '/' 'Male' '/' 'Noise' '/' m_nse_wav_fle{iFG}]);
    
    figure('Visible','off'); spectrogram(f_vce_wav, 128, [], [], fs_f_vce_wav, 'yaxis');
    z_lim = get(gca,'zlim');
    line('XData',[f_tme_hld{iFG,2} f_tme_hld{iFG,2}],'YData',[0 44100],'ZData',[z_lim(2) z_lim(2)],'LineWidth',0.4,'Color','w','LineWidth',2);
    line('XData',[f_tme_hld{iFG,3} f_tme_hld{iFG,3}],'YData',[0 44100],'ZData',[z_lim(2) z_lim(2)],'LineWidth',0.4,'Color','w','LineWidth',2);
    print(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Figures/Female/Voices' '/' f_vce_wav_fle{iFG}(1:end-4) '.png'],'-dpng')
    close all
    
    figure('Visible','off'); spectrogram(f_nse_wav, 128, [], [], fs_f_nse_wav, 'yaxis');
    z_lim = get(gca,'zlim');
    line('XData',[f_tme_hld{iFG,2} f_tme_hld{iFG,2}],'YData',[0 44100],'ZData',[z_lim(2) z_lim(2)],'LineWidth',0.4,'Color','w','LineWidth',2);
    line('XData',[f_tme_hld{iFG,3} f_tme_hld{iFG,3}],'YData',[0 44100],'ZData',[z_lim(2) z_lim(2)],'LineWidth',0.4,'Color','w','LineWidth',2);
    print(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Figures/Female/Noise' '/' f_nse_wav_fle{iFG}(1:end-4) '.png'],'-dpng')
    close all
    
%     figure('Visible','off'); spectrogram(m_vce_wav, 128, [], [], fs_m_vce_wav, 'yaxis');
%     z_lim = get(gca,'zlim');
%     line('XData',[m_tme_hld{iFG,2} m_tme_hld{iFG,2}],'YData',[0 44100],'ZData',[z_lim(2) z_lim(2)],'LineWidth',0.4,'Color','w','LineWidth',2);
%     line('XData',[m_tme_hld{iFG,3} m_tme_hld{iFG,3}],'YData',[0 44100],'ZData',[z_lim(2) z_lim(2)],'LineWidth',0.4,'Color','w','LineWidth',2);
%     print(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Figures/Male/Voices' '/' m_vce_wav_fle{iFG}(1:end-4) '.png'],'-dpng')
%     close all
%     
%     figure('Visible','off'); spectrogram(m_nse_wav, 128, [], [], fs_m_nse_wav, 'yaxis');
%     z_lim = get(gca,'zlim');
%     line('XData',[m_tme_hld{iFG,2} m_tme_hld{iFG,2}],'YData',[0 44100],'ZData',[z_lim(2) z_lim(2)],'LineWidth',0.4,'Color','w','LineWidth',2);
%     line('XData',[m_tme_hld{iFG,3} m_tme_hld{iFG,3}],'YData',[0 44100],'ZData',[z_lim(2) z_lim(2)],'LineWidth',0.4,'Color','w','LineWidth',2);
%     print(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Figures/Male/Noise' '/' m_nse_wav_fle{iFG}(1:end-4) '.png'],'-dpng')
%     close all
    
end

%% Make Grid Plots
ttt.con1 = con1;
ttt.con2 = con2;
ttt.con3 = con3;
ttt.vow1 = vow1;
ttt.vow2 = vow2;
ttt.vow3 = vow3;

f_tme_hld = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/female_phoneme_timing.csv'); [~,ind] = sort(f_tme_hld(:,1)); f_tme_hld = f_tme_hld(ind,:);

for iBL = 1:3;
    
    figure('Visible','off')
    
    cnt = 1;
    for iC = 1:numel(ttt.(['con' num2str(iBL)]))
        for iV = 1:numel(ttt.(['vow' num2str(iBL)]))
            
            stm = [ttt.(['con' num2str(iBL)]){iC} ttt.(['vow' num2str(iBL)]){iV}];
            
            wve_ind = find(strcmpi(f_vce_wav_fle,['F_' lower(stm) '.wav']));
            if isempty(wve_ind)
                wve_ind = find(strcmpi(f_vce_wav_fle,['F_' lower(stm(1:2)) '.wav']));
            end
            if isempty(wve_ind)
                wve_ind = find(strcmpi(f_vce_wav_fle,['F_' lower(stm(1:3)) '.wav']));
            end
            tme_ind = find(strcmpi([lower(stm)],f_tme_hld(:,1)));
            
            if iBL ~= 3; subplot(4,4,cnt); else subplot(4,4,cnt); end
            
            [f_vce_wav,fs_f_vce_wav] = wavread(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Stimuli/Female/Voices' '/' f_vce_wav_fle{wve_ind}]);
            
            spectrogram(f_vce_wav, 128, [], [], fs_f_vce_wav, 'yaxis');
            z_lim = get(gca,'zlim');
            ylabel([]); xlabel([]);
            line('XData',[f_tme_hld{tme_ind,2} f_tme_hld{tme_ind,2}],'YData',[0 44100],'ZData',[z_lim(2) z_lim(2)],'LineWidth',0.4,'Color','w','LineWidth',2);
            line('XData',[f_tme_hld{tme_ind,3} f_tme_hld{tme_ind,3}],'YData',[0 44100],'ZData',[z_lim(2) z_lim(2)],'LineWidth',0.4,'Color','w','LineWidth',2);
            
            title(stm)
            
            cnt = cnt+1;
        end
    end
    
    print(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/' '/' 'F_Block' num2str(iBL) '_Bigrams.png'],'-dpng','-r300');
    close all
    
end


for iBL = 1:3;
    
    figure('Visible','off')
    
    cnt = 1;
    for iC = 1:numel(ttt.(['con' num2str(iBL)]))
        for iV = 1:numel(ttt.(['vow' num2str(iBL)]))
            
            stm = [ttt.(['con' num2str(iBL)]){iC} ttt.(['vow' num2str(iBL)]){iV}];
            
            wve_ind = find(strcmpi(f_nse_wav_fle,['F_' lower(stm) '_matched_noise_TB12_Fc20-20000_nobroadband_log.wav']));
            if isempty(wve_ind)
                wve_ind = find(strcmpi(f_nse_wav_fle,['F_' lower(stm(1:2)) '_matched_noise_TB12_Fc20-20000_nobroadband_log.wav']));
            end
            if isempty(wve_ind)
                wve_ind = find(strcmpi(f_nse_wav_fle,['F_' lower(stm(1:3)) '_matched_noise_TB12_Fc20-20000_nobroadband_log.wav']));
            end
            tme_ind = find(strcmpi([lower(stm)],f_tme_hld(:,1)));
            
            if iBL ~= 3; subplot(4,4,cnt); else subplot(4,4,cnt); end
            
            [f_vce_wav,fs_f_vce_wav] = wavread(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/Stimuli/Female/Noise' '/' f_nse_wav_fle{wve_ind}]);
            
            spectrogram(f_vce_wav, 128, [], [], fs_f_vce_wav, 'yaxis');
            z_lim = get(gca,'zlim');
            ylabel([]); xlabel([]);
            line('XData',[f_tme_hld{tme_ind,2} f_tme_hld{tme_ind,2}],'YData',[0 44100],'ZData',[z_lim(2) z_lim(2)],'LineWidth',0.4,'Color','w','LineWidth',2);
            line('XData',[f_tme_hld{tme_ind,3} f_tme_hld{tme_ind,3}],'YData',[0 44100],'ZData',[z_lim(2) z_lim(2)],'LineWidth',0.4,'Color','w','LineWidth',2);
            
            title(stm)
            
            cnt = cnt+1;
        end
    end
    
    print(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/' '/' 'F_Block' num2str(iBL) '_Bigrams_noise.png'],'-dpng','-r300');
    close all
    
end

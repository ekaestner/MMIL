clear; clc;

tsk_dir = '/home/ekaestne/Desktop/Tasks/OrganizedTaskFolder/SL/SL_2_11_16';

low_pss = 13000;

Fs = 44100; % Hz
wrt_flg = 0; 
SCN_lab = 1; % whether to include word in SCN file name
out_str = '';
out_dir = '';
TB = 7;
Fc = [100 low_pss];

%% Fix lists
tsk_fle = dir([tsk_dir '/' 'lists' '/']); tsk_fle = tsk_fle(3:end);

for iL = 1:numel(tsk_fle)
    
    old_lst = mmil_readtext([tsk_dir '/' 'lists' '/' tsk_fle(iL).name]);
    
    for iR = 1:size(old_lst,1)
        if ~isempty(strfind(old_lst{iR,2},'_matched_noise_TB12_Fc20-20000_nobroadband_log.wav'))
            rmv_prt = strfind(old_lst{iR,2},'_matched_noise_TB12_Fc20-20000_nobroadband_log.wav');
            old_lst{iR,2} = old_lst{iR,2}(1:rmv_prt);
            old_lst{iR,2} = [old_lst{iR,2} 'nse.wav'];
        end 
        
    end
    
    cell2csv([tsk_dir '/' 'lists' '/' tsk_fle(iL).name],old_lst)
    
end

%% Make new noise
stm_fle = dir([tsk_dir '/' 'stimuli' '/']); stm_fle = stm_fle(3:end);

stm_cnt = 1;
for iS = 1:numel(stm_fle)
    
    if ~isempty(strfind(stm_fle(iS).name,'_matched_noise_TB12_Fc20-20000_nobroadband_log.wav'))
        
        delete([tsk_dir '/' 'stimuli' '/' stm_fle(iS).name])
        
    elseif isempty(strfind(stm_fle(iS).name,'_matched_noise_TB12_Fc20-20000_nobroadband_log.wav'))
     
        stm_hld{stm_cnt} = stm_fle(iS).name;
        stm_cnt          = stm_cnt + 1;       
        
        fem_fle = [tsk_dir '/' 'stimuli' '/' stm_fle(iS).name];
        [fem_nse_scn,word] = GenerateSCN_v2_remove2(fem_fle,Fc,TB,Fs,out_dir,[20 20],out_str,SCN_lab,iS,wrt_flg,1);
        fem_nse_scn        = fem_nse_scn*4;
        audiowrite([tsk_dir '/' 'stimuli' '/' stm_fle(iS).name(1:end-4) '_nse' stm_fle(iS).name(end-3:end)],fem_nse_scn,Fs);
        
    end
    
end

%% Make spectrograms to check noise
for iSP = 1:numel(stm_hld)
     if isempty(strfind(stm_hld{iSP},'_matched_noise_TB20_Fc50-5000.wav'))    
         
         [f_vce,f_vce_smp]         = audioread([tsk_dir '/' 'stimuli' '/' stm_hld{iSP}]);
         [f_vce_nse,f_vce_nse_smp] = audioread([tsk_dir '/' 'stimuli' '/' stm_hld{iSP}(1:end-4) '_nse' stm_hld{iSP}(end-3:end)]);
         
         figure('units','normalized','outerposition',[0 0 0.5 0.75],'Visible','off');
         subplot(1,2,1)
         spectrogram(f_vce, 128, [], [], f_vce_smp, 'yaxis'); colorbar;
         subplot(1,2,2)
         spectrogram(f_vce_nse, 128, [], [], f_vce_nse_smp, 'yaxis'); colorbar;
         
         print([tsk_dir '/' 'SL_Task_Info' '/' stm_hld{iSP} '.png'],'-dpng','-r250');
         
     end    
end

%% Fix noise
stm_fle = dir([tsk_dir '/' 'stimuli' '/']); stm_fle = stm_fle(3:end);

stm_cnt = 1;
for iS = 1:numel(stm_fle)
    
    if ~isempty(strfind(stm_fle(iS).name,'_nse.wav'))
        
        delete([tsk_dir '/' 'stimuli' '/' stm_fle(iS).name])
        
    elseif isempty(strfind(stm_fle(iS).name,'_matched_noise_TB12_Fc20-20000_nobroadband_log.wav'))
     
        stm_hld{stm_cnt} = stm_fle(iS).name;
        stm_cnt          = stm_cnt + 1;       
        
        fem_fle = [tsk_dir '/' 'stimuli' '/' stm_fle(iS).name];
        [fem_nse_scn,word] = GenerateSCN_v2_remove2(fem_fle,Fc,TB,Fs,out_dir,[20 20],out_str,SCN_lab,iS,wrt_flg,1);
        audiowrite([tsk_dir '/' 'stimuli' '/' stm_fle(iS).name(1:end-4) '_nse' stm_fle(iS).name(end-3:end)],fem_nse_scn,Fs);
        
    end
    
end







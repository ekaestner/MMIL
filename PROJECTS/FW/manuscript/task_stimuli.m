%% Repetition Stimuli
stm = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/task/fw_stimuli_2.csv');

rep_stm = find(cell2mat(stm(:,3))==4);
mean(diff(rep_stm))
range(diff(rep_stm))

%% Length of stimuli
wrd_stm = stm(find(cell2mat(stm(:,3))==3),4);
mean(cellfun(@numel,wrd_stm))
max(cellfun(@numel,wrd_stm))
min(cellfun(@numel,wrd_stm))

wrd_stm = stm(find(cell2mat(stm(:,3))==4),4);
mean(cellfun(@numel,wrd_stm))
max(cellfun(@numel,wrd_stm))
min(cellfun(@numel,wrd_stm))

wrd_stm = stm(find(cell2mat(stm(:,3))==5),4);
mean(cellfun(@numel,wrd_stm))
max(cellfun(@numel,wrd_stm))
min(cellfun(@numel,wrd_stm))

wrd_stm = stm(find(cell2mat(stm(:,3))==6),4);
mean(cellfun(@numel,wrd_stm))
max(cellfun(@numel,wrd_stm))
min(cellfun(@numel,wrd_stm))

wrd_stm = stm(find(cell2mat(stm(:,3))==7),4);
mean(cellfun(@numel,wrd_stm))
max(cellfun(@numel,wrd_stm))
min(cellfun(@numel,wrd_stm))

%% Stimulus Characteristics
stm_loc = stm(find(cell2mat(stm(:,3))~=7),:);
stm_chr = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/task/fw_2.csv');

stm_ort = cell2mat(stm_chr(find(cell2mat(stm_loc(:,3))==3),4));
mean(stm_ort)
min(stm_ort)
max(stm_ort)

stm_frq = cell2mat(stm_chr(find(cell2mat(stm_loc(:,3))==3),3));
nanmean(stm_frq)
min(stm_frq)
max(stm_frq)

stm_big = cell2mat(stm_chr(find(cell2mat(stm_loc(:,3))==3),5));
mean(stm_big)
min(stm_big)
max(stm_big)

%%
ttt = load('//space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/machinelearning/WRD_FFN/regions/middle-precentral_L_acc.mat');

ttt.reg.acc




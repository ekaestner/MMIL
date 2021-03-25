%% Make lists for ML

Con = {'' '' '' '' '' '' '' '' '' ''};
Vow = {'' '' '' '' '' '' '' '' '' ''};

Stim_Num = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];

stl = length(A_Words)*(length(A_Words)-1);
stim_list = cell(stl*4,4);

% Create Audio Versions of A_Words
for ic = 1:length(A_Words_norm)
    A_Words{ic} = ['F_' A_Words_norm{ic} '.wav'];
end

% Create Auditory Noise Versions of A_Words
for ic = 1:length(A_Words_norm)
    A_noise_Words{ic} = ['F_' A_Words_norm{ic} '_matched_noise_TB20_Fc50-5000.wav'];
end

% matched_words
for iword = 0:length(A_Words)-1
    
    stim_list(iword*12+1:(iword+1)*12,1) = {A_Words{iword+1}};
    stim_list(iword*12+1:(iword+1)*12,2) = {1};
    stim_list(iword*12+1:(iword+1)*12,3) = {Stim_Num(iword+1)};
    stim_list(iword*12+1:(iword+1)*12,4) = {Stim_Num(iword+1)};
    
end

% mismatch_words
for imis = 0:length(A_Words)-1
    
    stim_list(stl+imis*12+1:stl+(imis+1)*12,1) = {A_Words{imis+1}};
    stim_list(stl+imis*12+1:stl+(imis+1)*12,2) = {2};
    stim_list(stl+imis*12+1:stl+(imis+1)*12,3) = {Stim_Num(imis+1)};
    
    counter = stl+imis*12+1;
    
    for iv = 1:length(A_Words)
        if iv ~= imis+1
            stim_list(counter,4) = {Stim_Num(iv)};
            counter = counter + 1;
        end
    end
end

% Auditory Noise - Visual Word
for imis = 0:length(A_Words)-1
    
    stim_list(stl*2+imis*12+1:stl*2+(imis+1)*12,1) = {A_noise_Words{imis+1}};
    stim_list(stl*2+imis*12+1:stl*2+(imis+1)*12,2) = {3};
    stim_list(stl*2+imis*12+1:stl*2+(imis+1)*12,3) = {Stim_Num(imis+1)};
    
    counter = stl*2+imis*12+1;
    
    for iv = 1:length(A_Words)
        if iv ~= imis+1
            stim_list(counter,4) = {Stim_Num(iv)};
            counter = counter + 1;
        end
    end
end

% Auditory Word - Visual Noise
for imis = 0:length(A_Words)-1
    
    stim_list(stl*3+imis*12+1:stl*3+(imis+1)*12,1) = {A_Words{imis+1}};
    stim_list(stl*3+imis*12+1:stl*3+(imis+1)*12,2) = {4};
    stim_list(stl*3+imis*12+1:stl*3+(imis+1)*12,3) = {Stim_Num(imis+1)};
    
    counter = stl*3+imis*12+1;
    
    for iv = 1:length(A_Words)
        if iv ~= imis+1
            stim_list(counter,4) = {Stim_Num(iv)};
            counter = counter + 1;
        end
    end
end

% Save whole list
cell2csv('/home/ekaestne/Desktop/cvcvc/SL/stimlist.csv',stim_list);

% Randomize list
stim_list_randomized = ejk_createlist(stim_list,[3 1 1],3);
cell2csv('/home/ekaestne/Desktop/cvcvc/SL/stimlist_randomized.csv',stim_list_randomized);

% Make Audio List
Alist = stim_list_randomized(:,[1 2 3]);

% Make Visual List
Vlist(:,2:3) = stim_list_randomized(:,[2 4]);
for ipop = 1:length(Vlist)
    Vlist{ipop,1} = V_Words{Vlist{ipop,3}};
end
    
% Save List in chunks
ejk_presentationlist(Alist,'/home/ekaestne/Desktop/cvcvc/SL/','AList',6)
ejk_presentationlist(Vlist,'/home/ekaestne/Desktop/cvcvc/SL/','VList',6)



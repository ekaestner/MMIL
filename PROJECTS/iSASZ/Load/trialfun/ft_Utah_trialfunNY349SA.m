%% Function to read events & create file to epoch NYU iEEG data from a edf.mat format
% Function to epoch NYU iEEG data. 
% dat has 8 channels:each channel has the same trigger patten as the others
% Created by Jane Deng (4-13-15) q5deng@ucsd.edu

function [trl, evt] = ft_Utah_trialfunNY349SA(data)

addpath(genpath('/home/qdeng/qdeng'));
addpath(genpath('/space/syn02/1/data/MMILDB/BACKUP/md7/1/pubsw/packages/MMPS/MMPS_236'));

data = load('/space/mdeh4/1/halgdev/projects/mmilanguage/TriggerChannel/NY349_SA_edf.mat');
data_orig_name = cell2mat(fieldnames(data));
dat = getfield(data,data_orig_name);

thresh = 0.2 ;

Fs = 512;
pre = 1.5*Fs;
pos = 2.5*Fs;
%% Find the timings
tmpdat = dat.trial{1}(8,:);   tmpdat([1:14e4])=0;  tmpdat([3.7e5:2.38e6])=0;   
   
ind     = crossing(tmpdat,[],thresh);  
ind     = ind(1:2:length(ind));         % only keep left edge of pulse only the odd number
tmpdat2 = tmpdat(ind+2);                % select the peak of each pulses

trig_pos   =  [1 find(diff(ind)*1000/Fs > 100)+1];
samp   = ind(trig_pos);
% if sum(samp-pre < 1) || sum(samp+pos > numel(dat.time{1}));%delete the remote trial, and
%    trig_pos(samp+pos > numel(dat.time{1})) =[];
%    trig_pos(samp-pre < 1) =[];
%    samp(samp+pos > numel(dat.time{1})) =[];
%    samp(samp-pre < 1) =[];
% end

Nevents = numel(samp);
sel = zeros(1,numel(ind));
sel(trig_pos) = 1;
r = round(tmpdat2);
trigs = cell(1,Nevents);

%ft
for i = 1:Nevents;
    trigs{i} = r([trig_pos(i)]);
    trig_num(i) = numel(trigs{i});
end
% trigs{Nevents} = r(trig_pos(Nevents)+1:end);
% trig_num(Nevents) = numel(trigs{Nevents});

% if sum(ismember(trig_num,0)) >= 1 
%     ismember(trig_num,0);
% end
% 
%  z_ind = ismember(trig_num,0);
% % trigs(z_ind) = 1;              % throw away 0 member in trigger information
% % samp(z_ind) = [];               % throw away 0 member in sample
%  trig_num(z_ind) = 1;           % throw away 0 member in trigger number information 

%% Define Events
evt = zeros(1,numel(trig_num));

trigdata = mmil_readtext('/home/qdeng/qdeng/dataset/SAi-actual.csv');
trig =trigdata(:,3);
evt_val = cell2num(trig);

%% Adjust time to make it in the correct sampling frequency

% this one has been downsampled already.no need to run this step
%     dsfact = 3e4/500;
%     samp  = round(samp/dsfact);  
%     pre = round(pre/dsfact);
%     pos = round(pos/dsfact);


%% Make new trl file
mtx = zeros(length(evt),4);
for i=1:length(evt)
    mtx(i,1) = samp(i)-pre;
    mtx(i,2) = samp(i)+pos;
    mtx(i,3) = -pre;
    mtx(i,4) = evt(i);
end
tri_end = size(mtx,1);
tri_blc = find(diff(mtx(:,1))>3000);

tri_blc_1 = [1:1:tri_blc(1)];
tri_blc_2 = [(tri_blc(1)+1):1:tri_blc(2)];
tri_blc_3 = [(tri_blc(2)+1):1:tri_blc(3)];
tri_blc_4 = [(tri_blc(3)+1):1:tri_blc(4)];
tri_blc_5 = [(tri_blc(4)+1):1:tri_end];


%  for block 1
for i = 1:length(tri_blc_1)
    sen1(i) = abs(mtx(2,1)-mtx(tri_blc_1(i),1));
end
sen2 = zeros(1,length(sen1)-1);
for i = 3:(length(sen1)-2)
    if sen1(i+1)-sen1(i) > 1100
       sen2(i) =sen1(i);
       sen2(i+1)= sen1(i+1);
    elseif (sen1(i+1)-sen1(i) < 1000) && sen2(i) ~= 0 && (sen1(i+2)-sen1(i) > 1000)
            sen2(i+2) = sen1(i+2);      
    end  
   if sen1(end)-sen1(end-1) > 1100
      sen2(end) = sen1(end);
   end
end
idd = find(~ismember(sen2,0));
if mtx(length(tri_blc_1),1)-mtx(length(sen2),1) >1000
   idd(end+1) = length(tri_blc_1);
end
seq1= sort([tri_blc_1(2) tri_blc_1(idd)]);


% for block 2
for i = 1:length(tri_blc_2)
    sen3(i) = abs(mtx(tri_blc_2(2),1)-mtx(tri_blc_2(i),1));
end
sen4 = zeros(1,length(sen3)-1);
for i = 2:(length(sen3)-2)
    if sen3(i+1)-sen3(i) > 1100
       sen4(i) =sen3(i);
       sen4(i+1)= sen3(i+1);
    elseif (sen3(i+1)-sen3(i) < 1000) && sen4(i) ~= 0 && (sen3(i+2)-sen3(i) > 1000)
            sen4(i+2) = sen3(i+2);      
    end  
   if sen3(end)-sen3(end-1) > 1100
      sen4(end) = sen3(end);
   end
end
idd = find(~ismember(sen4,0));
if mtx(tri_blc_2(end),1)-mtx(tri_blc_2(end)-1,1) >1000
   idd(end+1) = idd(end)+1;
end
seq2= sort([tri_blc_2(2) tri_blc_2(idd) 420 422 424 426]);


% for block 3

for i = 1:length(tri_blc_3)
    sen5(i) = abs(mtx(tri_blc_3(2),1)-mtx(tri_blc_3(i),1));
end
sen6 = zeros(1,length(sen5)-1);
for i = 2:(length(sen5)-2)
    if sen5(i) == 0
       sen6(i) = 1;
    end
    if sen5(i+1)-sen5(i) > 1000
       sen6(i) =sen5(i);
       sen6(i+1)= sen5(i+1);
    elseif (sen5(i+1)-sen5(i) < 1000) && sen6(i) ~= 0 && (sen5(i+2)-sen5(i) > 1000)
            sen6(i+2) = sen5(i+2);      
    end  
    if sen5(end)-sen5(end-1) > 1000
       sen6(end) = sen5(end);
    end
end
idd = find(~ismember(sen6,0));
if mtx(tri_blc_3(end),1)-mtx(tri_blc_3(end)-1,1) >1000
   idd(end+1) = idd(end)+1;
end
seq3= sort([tri_blc_3(idd)]); 
rid = find(ismember(seq3,512)); seq3(rid)=[];


% for block 4

for i = 1:length(tri_blc_4)
    sen7(i) = abs(mtx(tri_blc_4(2),1)-mtx(tri_blc_4(i),1));
end
sen8 = zeros(1,length(sen7)-1);
for i = 2:(length(sen7)-2)
    if sen7(i+1)-sen7(i) > 1000
       sen8(i) =sen7(i);
       sen8(i+1)= sen7(i+1);
    elseif (sen7(i+1)-sen7(i) < 1000) && sen8(i) ~= 0 && (sen7(i+2)-sen7(i) > 1000)
            sen8(i+2) = sen7(i+2);      
    end  
    if sen7(end)-sen7(end-1) > 1000
       sen8(end) = sen7(end);
    end
end
idd = find(~ismember(sen8,0));
if mtx(tri_blc_4(end),1)-mtx(tri_blc_4(end)-1,1) >1000
   idd(end+1) = idd(end)+1;
end
seq4= sort([tri_blc_4(2) tri_blc_4(idd)]);


% for block 5

for i = 1:length(tri_blc_5)
    sen9(i) = abs(mtx(tri_blc_5(2),1)-mtx(tri_blc_5(i),1));
end
sen10 = zeros(1,length(sen9)-1);
for i = 2:(length(sen9)-2)
    if sen9(i+1)-sen9(i) > 1000
       sen10(i) =sen9(i);
       sen10(i+1)= sen9(i+1);
    elseif (sen9(i+1)-sen9(i) < 1000) && sen10(i) ~= 0 && (sen9(i+2)-sen9(i) > 1000)
            sen10(i+2) = sen9(i+2);      
    end  
    if sen9(end)-sen9(end-1) > 1000
       sen10(end) = sen9(end);
    end
end
idd = find(~ismember(sen10,0));
if mtx(tri_blc_5(end),1)-mtx(tri_blc_5(end)-1,1) >1000
   idd(end+1) = idd(end)+1;
end
seq5= sort([tri_blc_5(2) tri_blc_5(idd)]);


%% make trigger matrix
id = sort([seq1 seq2 seq3 seq4 seq5]);
trl=mtx(id,:);
trl(:,4) = evt_val(1:744);

rmpath(genpath('/home/qdeng/qdeng'));
rmpath(genpath('/space/syn02/1/data/MMILDB/BACKUP/md7/1/pubsw/packages/MMPS/MMPS_236'));

end



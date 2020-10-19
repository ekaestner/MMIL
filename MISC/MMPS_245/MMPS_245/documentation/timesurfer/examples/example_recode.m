% test_recode

matfile = 'data/proc_events_2.mat';
load(matfile);

%recode_rules = {'2->7=101','15->8=102','8->15->6=103'};
%recode_rules = {'[2]->200:1100[7,18]=101','15->8=102','8->15->6=103'};
%recode_rules = {'[2,15]->7=101'};
%recode_rules = {'2->7->15=101'};
%recode_rules = {'2->200:1100[15]=101'};
%recode_rules = {'15->8=102','102->8=110'};
recode_rules = {'2->15=101','2->7=102'};

[mod_evnts,new_evnts] = ts_recode_events(evnts,recode_rules,hdr.sfreq);



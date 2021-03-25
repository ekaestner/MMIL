function varargout = ts_iEEG_findevents (eegfile)
% Usage:  ts_iEEG_findevents ('eeg file name')
%
% prints out the list of event types from an eeg file
%
% Created on ??? by Rajan Patel

data = pop_loadeeg(eegfile);

for i=1:data.trials
        conds(i)=data.event(i).epochtype;
end

types=unique(conds);
fprintf('\nEvent types in %s:\n',eegfile)
for i=1:length(types)
    fprintf('%s\n',num2str(types(i)));
end

if nargout > 0
  varargout{1} = types;
end
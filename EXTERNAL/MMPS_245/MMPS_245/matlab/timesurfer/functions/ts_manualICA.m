function epoch_data = ts_manualICA(epoch_data,varargin)
%function epoch_data = ts_manualICA (epoch_data,'maxsteps',30,...)
%
% Usage:
%  epoch_data = ts_manualICA (epoch_data,'key',value,...);
%
% Required Input:
%  epoch_data  - epoch data as read by ts_process_fif_data
%
% Optional Input:
%  'maxsteps':  maximum number of steps for ica
%     {default = 20}
%  'eps': minimum amplitude
%     {default = 1e-70}
%  'plottype: type of plot to view 'alltrials' or 'activations'
%    'activations': shows from 1 to 'ntrial' trials
%                   of component activations along x-axis
%    'alltrials': shows all trials (y-axis) as a color plot
%                 with time on x-axis and activation size on c-axis
%     {default = 'activations'}
%  'ntrial': the number of trials to display on screen for IC
%            selection; only applies to 'activations' plot
%     {default = 5}
%  'start_trial': first trial to display on screen for IC
%                 selection; only applies to 'activations' plot
%     {default = 1}
%  'ncomponents': the number of components to display on screen for IC
%                  selection; use [] to display all
%     {default = 80}
%  'compperfig': number of components to display per figure
%     {default = 40}
%  'event_codes'
%       OR
%  'conditions': the event codes or conditions to inspect
%     {default = 'all'}
%  'chantype': cell array of channel types to process
%     each set will be processed individually
%  'outfile': if called without an output data saved to this file
%
% Output:
%
%   epoch_data: epoch data without the selected components
%
%  If no output specfied script will save epoch_data to outfile.
%
% Created:   03/18/08 by Andrei Irimia
% Modified:  06/26/08 by Matt Leonard
% Modidied:  09/01/11 by Burke Rosen -added start_trial
% Last Mod:  08/11/15 by Don Hagler
%

%% Check Inputs

if nargin < 1, help(mfilename); end

parms = mmil_args2parms(varargin,...
                        { ...
                         'ntrial',5,[],...
                         'start_trial',1,[],...
                         'plottype','activations',{'alltrials','activations'},...
                         'ncomponents',80,[],...
                         'compperfig',40,[],...
                         'maxsteps',20,[],...
                         'eps',1e-70,[],...
                         'event_codes',[],[],...
                         'conditions',[],[],...
                         'chantype', 'all', {'all', 'mag' 'grad1' 'grad2' 'eeg', 'other', 'grad', 'meg'},...
                         'outfile',[],[],...
                         'verbose',true,sort([false true]),...
                         'logfile',[],[],...
                         'logfid',1,[],...
                         'rescale',true,sort([false true]),...
                         'rootoutdir',pwd,[],...
                         'prefix','proc',[],...
                         'ICA_saveout_flag',1,{0,1}...
                        },...
                        false);

epoch_data = ts_checkdata_header(epoch_data);

% store pre-ica epoch data in order to scale post-ICA data to original means and variances
pre_ICA_epoch_data = epoch_data;

errors = ts_checkdata(epoch_data);
if ~isempty(errors)
    mmil_error(parms,'Errors in provided data structure: %s.',errors);
end                   

if (~isempty(parms.event_codes) && ~isempty(parms.conditions))
    mmil_error(parms, 'Cannot specify BOTH event_codes AND conditions.');

    %specified none
elseif (isempty(parms.event_codes) && isempty(parms.conditions))
    parms.event_codes = { epoch_data.epochs.event_code };
    parms.conditions  = num2cell(1:length(epoch_data.epochs));

    % specified conditions
elseif (isempty(parms.event_codes))
    if (min(parms.conditions) < 1 || max(parms.conditions) > length(epoch_data.epochs))
        mmil_error(parms, 'Conditions are out of range; nConditions=%d', length(epoch_data.epochs));
    else
        parms.event_codes = { epoch_data.epochs(parms.conditions).event_code };
    end;

    % specified event_codes
else
    if (~isempty(setdiff([parms.event_codes{:}], [epoch_data.epochs.event_code])))
        mmil_error(parms, 'Event code doesn''t exist in epoch data: %d.', ...
            min(setdiff([parms.event_codes{:}], [epoch_data.epochs.event_code])));
    else
        [a,parms.conditions]=intersect([epoch_data.epochs.event_code], [parms.event_codes{:}]);
        parms.conditions    = num2cell(parms.conditions);
    end;
end;

% Make sure both conditions and event_codes are cell arrays
if (~iscell(parms.event_codes)), parms.event_codes = num2cell(parms.event_codes); end;
if (~iscell(parms.conditions)),  parms.conditions  = num2cell(parms.conditions);  end;

if ~isempty(parms.chantype)
    if ~iscell(parms.chantype), parms.chantype = {parms.chantype}; end
    selchantypes = {};
    for c = 1:length(parms.chantype)
        switch parms.chantype{c}
            case 'grad1'
                selchantypes{end+1} = 'grad1';
            case 'grad2'
                selchantypes{end+1} = 'grad2';
            case 'grad'
                selchantypes{end+1} = 'grad';
            case 'all'
                selchantypes = {'grad','mag','eeg'};
            case 'mag'
                selchantypes{end+1} = 'mag';
            case 'meg'
                selchantypes{end+1} = 'mag';
                selchantypes{end+1} = 'grad';
            case 'other'
                selchantypes{end+1} = 'other';
            case 'eeg'
                selchantypes{end+1} = 'eeg';
        end
    end
else
    selchantypes = {'grad','mag','eeg'};
end

%%  Run ICA     

% % find bad channels
bchans  = find([epoch_data.sensor_info.badchan]);

for k = 1:length(parms.conditions)
    [numchan,samples,numtrials] = size(epoch_data.epochs(parms.conditions{k}).data);
    data_in = reshape(epoch_data.epochs(parms.conditions{k}).data,numchan,samples*numtrials);
    for chtypes = selchantypes % go through each channel type separately
        chans = [];
        switch chtypes{1}
            case {'grad1','grad2','mag' 'eeg', 'other'}
                chans = find(strcmp(chtypes{1},{epoch_data.sensor_info.typestring}));
            case {'grad'}
                chans = find(strncmp(chtypes{1},{epoch_data.sensor_info.typestring},...
                    length(chtypes{1})));
        end;
        chans = setdiff(chans,bchans); % ignore bad channels
        if ~isempty(chans)
            mmil_logstr(parms,'Processing %s channels for event %d.',chtypes{1},epoch_data.epochs(parms.conditions{k}).event_code);
            dataset = data_in(chans,:);
            cnt = 0; bad = [];
            for m = 1:size(dataset,1)
                if abs(max(dataset(m,:))) < parms.eps && abs(min(dataset(m,:))) < parms.eps
                    cnt = cnt + 1;
                    bad(cnt) = m;
                else
                    if median(dataset(m,:)) ~= 0, dataset(m,:) = dataset(m,:)./abs(median(dataset(m,:))); end
                end
            end
            for d = cnt:-1:1, dataset(bad(d),:) = []; end
            chans(bad) = [];
            
            %% if ICA doesn't run
            weird = 0;
            if weird
                upper = 6;
                for d = 1:upper
                    len = size(dataset,2)/upper;
                    dta = dataset(:,(d-1)*len + 1:d*len);
                    [weights, sphere] = runica(dta,'maxsteps',parms.maxsteps,'verbose','on');
                    IC(d,:,:) = weights*sphere*dta;
                end
                IC = permute(IC,[2 1 3]);
                IC = reshape(IC,size(IC,1),size(IC,2)*size(IC,3));
            else % for normal ICA
                [weights, sphere] = runica(dataset,'maxsteps',parms.maxsteps,'verbose','on');
                IC = weights*sphere*dataset;
            end
            % Sort by latency as well ???
            % [windex,maxvar,maxframe,maxepoch,maxmap] = compsort(dataset,weights,sphere,mean(dataset')');
            % weights = weights(windex,:);
            
            if ~isempty(parms.ncomponents)
                parms.ncomponents = min(parms.ncomponents, size(IC,1));
            else
                parms.ncomponents = size(IC,1);
            end
            numplots = ceil(parms.ncomponents/parms.compperfig);
            % Display ICs
            t  = 0:1/samples:parms.ntrial-1/samples;
            for p = 1:numplots
                figure('CloseRequestFcn',@my_closefcn,'Units','normalized','Position',[0 0 1 1]);
                for m = (parms.compperfig*(p-1))+1:min((parms.compperfig*(p-1))+parms.compperfig,size(IC,1))
                    complegend(m).handle    = subplot(10,ceil(parms.compperfig/10),m-(parms.compperfig*(p-1)));
                    switch parms.plottype
                        case 'activations', plot(t,IC(m,parms.start_trial*samples:((parms.start_trial+parms.ntrial)*samples)-1));
                        case 'alltrials', imagesc(permute(reshape(IC(m,:),samples,numtrials),[2 1]));
                    end
                    axis tight;
                    set(gca,'xticklabel',num2str(parms.start_trial+str2num(get(gca,'xticklabel'))));
                    title([' IC #' int2str(m)],'Color','k');
                    drawnow;
                end
            end
            % Select Bad Components
            selection = 'No';
            R        = [];
            tempcomp = [];
            while strcmpi(selection,'No');
                fprintf('%s: Toggle the components you wish to remove with the mouse.\n',mfilename);
                fprintf('%s: Press any key to end.\n',mfilename);
                b = waitforbuttonpress;
                while b == 0 % continue until keyboard button is pressed
                    selcomph   = get(gcf,'CurrentAxes');  % Axes handle
                    selpos     = get(gcf,'CurrentPoint'); % [x y]
                    selcomppos = get(selcomph,'Position');% actual position of subplot [left bottom width height]
                    if selpos(1) >= selcomppos(1) && selpos(1) <= (selcomppos(1)+selcomppos(3)) && ...  % make sure the selection
                            selpos(2) >= selcomppos(2) && selpos(2) <= (selcomppos(2)+selcomppos(4))         % is in range -> this can happen in the final figure when there are only a couple subplots
                        tempcomp   = find([complegend.handle] == selcomph); % find matching component
                        if find(R == tempcomp)                              % toggle selected component
                            R(find(R == tempcomp)) = [];                                     % remove from list
                            title(selcomph,[' IC #' int2str(tempcomp)],'Color','k');         % restore to normal
                        else
                            R(end+1) = tempcomp;                                             % add to list
                            title(selcomph,['IC #' int2str(tempcomp) '\surd'],'Color','r'); % color red with a check mark
                        end
                        R = sort(R);                                                       % sort
                        tempcomp = [];
                    end
                    b = waitforbuttonpress;
                end
                if ~isempty(R)
                    fprintf('%s: You have currently selected the following components:\n',mfilename);
                    fprintf('%s: %s\n',mfilename,regexprep(num2str(R),'\s*',', '));
                else
                    fprintf('%s: You have not selected any components to remove.\n',mfilename);
                end
                selection = questdlg('Are you happy with your selections?',...
                    'Done?',...
                    'Yes','No','No');
            end
            close all force;
            mmil_logstr(parms,'Removing %d selected component(s) from the %s channels.',length(R),chtypes{1});
            reject_data = icaproj(dataset,weights*sphere,setdiff(1:size(IC,1),R),mean(dataset')');
            epoch_data.epochs(parms.conditions{k}).data(chans,:,:) = reshape(reject_data,size(reject_data,1),samples,numtrials);
            clear weights sphere dataset reject_data
        else
            mmil_logstr(parms,'No good %s channels found.',chtypes{1});
        end
    end
    data_in = [];
end

% if ~exist(sprintf('%s/matfiles',parms.rootoutdir))
%   mkdir(sprintf('%s/matfiles',parms.rootoutdir));
% end


% rescale data after ICA to recover mean and variance
% if ~exist(sprintf('%s/matfiles/%s_epoch_data.mat',parms.rootoutdir,parms.prefix));
%     rawfile = sprintf('%s/matfiles/%s_epoch_data_1.mat',parms.rootoutdir,parms.prefix); % for testing purposes
% else
%     rawfile = sprintf('%s/matfiles/%s_epoch_data.mat',parms.rootoutdir,parms.prefix);
% end
% icafile = sprintf('%s/matfiles/%s_epoch_data_ICA.mat',parms.rootoutdir,parms.prefix);
outfile = sprintf('%s/matfiles/%s_epoch_data_ICA.mat',parms.rootoutdir,parms.prefix);
% epoch_data = ts_ICAfit_manual(rawfile,icafile,outfile,'prefix',parms.prefix,'rootoutdir',parms.rootoutdir);
if parms.ICA_saveout_flag
    save(sprintf('%s/matfiles/%s_epoch_data_ICA.mat',parms.rootoutdir,parms.prefix),'epoch_data','-v7.3');
    epoch_data = ts_ICAfit_manual(pre_ICA_epoch_data,epoch_data,'outfile',outfile,'prefix',parms.prefix,'rootoutdir',parms.rootoutdir);
else
    epoch_data = ts_ICAfit_manual(pre_ICA_epoch_data,epoch_data);
end

fprintf('%s: manual ICA finished.\n',mfilename);

function my_closefcn(src,evnt)
% Do nothing, closing a figure
% will crash waitforbuttonpress



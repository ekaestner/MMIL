function outdata = ts_equate_conds(indata,varargin)
%Purpose: Removes trials at random until all conditions have an equal
%         number of trials. Unless otherwise specified (see below) the
%         target number of trials is equal the number of trials in the
%         condition with the fewest trials. 
%
%Example: epoch_data = ts_equate_conds(epoch_data,'tolerance',5)
%
%Inputs: timesurfer epoch data or timefreq data structure.
%Outputs: timesurfer structure of the same types with equated trials.
%Parameters: default : options :  description
%   maxtrialcount: min(indata.(datafield).num_trials) : [] : maximum nuber of trials per condition 
%   tolerance: 0 : [] : number of trials by which a codition may excede the maxtrialcount
%   
%Created by BQR 06/11/12

parms = mmil_args2parms( varargin, ...
    {'tolerance',0,[],...
     'maxtrialcount',[],[],...
    },false );
outdata = indata; clear indata;
if isfield(outdata,'epochs')
    if isempty(parms.maxtrialcount)
        parms.maxtrialcount = min([outdata.epochs.num_trials]);
    end
    for ieve = 1:length(outdata.epochs)
        while outdata.epochs(ieve).num_trials > parms.maxtrialcount+parms.tolerance;
          equat_idx = randi(outdata.epochs(ieve).num_trials);
          outdata.epochs(ieve).num_trials = outdata.epochs(ieve).num_trials - 1;
          fn = fieldnames(outdata.epochs(ieve).trial_info);
          for ifn = 1:length(fn)
             outdata.epochs(ieve).trial_info.(fn{ifn})(equat_idx) = [];
          end
          outdata.epochs(ieve).data(:,:,equat_idx) = [];
        end 
    end    
elseif isfield(outdata,'timefreq')
    if isempty(parms.maxtrialcount)
        parms.maxtrialcount = min([outdata.timefreq.num_trials]);
    end
    for ieve = 1:length(outdata.timefreq)
        while outdata.timefreq(ieve).num_trials > parms.maxtrialcount+parms.tolerance;
          equat_idx = randi(outdata.timefreq(ieve).num_trials);
          outdata.timefreq(ieve).num_trials = outdata.timefreq(ieve).num_trials - 1;
          fn = fieldnames(outdata.timefreq(ieve).trial_info);
          for ifn = 1:length(fn)
             outdata.timefreq(ieve).trial_info.(fn{ifn})(equat_idx) = [];
          end
          if isfield(outdata.timefreq,'power')
            outdata.timefreq(ieve).power(:,:,:,equat_idx) = [];
          end
          if isfield(outdata.timefreq,'cmplx')
            outdata.timefreq(ieve).cmplx(:,:,:,equat_idx) = [];
          end
        end
    end  
end
end 
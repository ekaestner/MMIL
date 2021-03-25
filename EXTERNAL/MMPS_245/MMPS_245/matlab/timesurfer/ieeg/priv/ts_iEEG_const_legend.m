function [figure_legend] = ts_iEEG_const_legend (num_legends, data_info, data_procs)

%  Creates a figure legends based on information provided by
%  ts_iEEG_decon_fname.
%
%  Created 08/28/2007 - Rajan H Patel
% 
%  See also: ts_iEEG_decon_fname

figure_legend = {};

for i = 1:num_legends
    bipolar     = '';
    avg         = '';
    tel         = '';
    lam         = '';
    low_pass    = '';
    high_pass   = '';
    notch       = '';
    fil_type    = '';
    fil_db      = '';
    bc          = '';
    rej         = '';
    dt          = '';
    ln          = '';
    for j = 1:size(data_info,2)
        switch char(data_info(i,j))
            case 'bipolar'
                bipolar = sprintf('Bipolar\n');
            case'avg'
                avg     = sprintf('Average\n');
            case 'tel'
                tel     = sprintf('Telemetry\n');
            case 'lam'
                lam     = sprintf('Laminar\n');
        end
    end
    if ~isempty(data_procs{1})
     for j = 1:size(data_procs,2)
       if ~isempty(data_procs{i,j})
        if findstr(char(data_procs(i,j)),'lp'), low_pass   = sprintf('Low Pass Filter: %sHz\n',strtok(char(data_procs(i,j)),'lp')); end;
        if findstr(char(data_procs(i,j)),'hp'), high_pass  = sprintf('High Pass Filter: %sHz\n',strtok(char(data_procs(i,j)),'hp'));end;
        if findstr(char(data_procs(i,j)),'fil'),fil_type   = sprintf('Filter Type: %s\n',strtok(char(data_procs(i,j)),'fil'));      end;
        if findstr(char(data_procs(i,j)),'db'), fil_db     = sprintf('Filter db: %sdb\n',strtok(char(data_procs(i,j)),'db'));       end;
        if findstr(char(data_procs(i,j)),'nf'), notch      = sprintf('Notch Filter\n');                                             end;
        if findstr(char(data_procs(i,j)),'bc'), bc         = sprintf('%s\n','Baseline Corrected');                                  end; 
        if findstr(char(data_procs(i,j)),'dt'), dt         = sprintf('%s\n','Detrend On');                                          end;
        if findstr(char(data_procs(i,j)),'rej'),rej        = sprintf('%s\n','Manual Rejection');                                    end; 
        if findstr(char(data_procs(i,j)),'ln'), ln         = sprintf('%s\n','Line Filter');                                         end;
       end
     end
    end
    figure_legend(i) = {sprintf('%s%s%s%s%s%s%s%s%s%s%s%s%s',bipolar,avg,tel,lam,low_pass,high_pass,fil_type,fil_db,notch,bc,dt,rej,ln)};
end
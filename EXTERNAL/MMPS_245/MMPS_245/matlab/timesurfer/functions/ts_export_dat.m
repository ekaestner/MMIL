function ts_export_dat(indata,outfname)
% ts_export_dat(indata,outfname)
% exports continuous data in timesurfer format to .dat ASCII file readable
% by neuroscan. 
% Required paramters:
% indata: Timesurfer datastruct containing continous data
% outfname: fullfile path to output file
%
% Created 02/19/13 BQR
% TODO: detect and accept average and epoch data inputs

fid = fopen(outfname,'w');
%% header
fprintf(fid, '[Subject]\txxxxx, xxxxx\r\n');
fprintf(fid, '[Date]\r\n');
fprintf(fid, '[Time]\r\n');
fprintf(fid, '[Channels]\t%i\r\n',indata.num_sensors);
fprintf(fid, '[Rate]\t %i\r\n',indata.sfreq);
fprintf(fid, '[Type]\tContinuous\r\n');
fprintf(fid, '[Rows]\tPoints\r\n');
fprintf(fid, '[Electrode Labels]\r\n');
for ichan = 1:indata.num_sensors
    fprintf(fid,'[%s]',indata.sensor_info(ichan).label);
    if ichan<indata.num_sensors
        fprintf(fid,'\t');
    end
end;fprintf(fid,'\r\n');
fprintf(fid, '[Electrode XUnits]\r\n');
for ichan = 1:indata.num_sensors
    fprintf(fid,'[Default]');
    if ichan<indata.num_sensors
        fprintf(fid,'\t');
    end
end;fprintf(fid,'\r\n');
fprintf(fid, '[Electrode YUnits]\r\n');
for ichan = 1:indata.num_sensors
    fprintf(fid,'[uV]');
    if ichan<indata.num_sensors
        fprintf(fid,'\t');
    end
end;fprintf(fid,'\r\n');
fprintf(fid, '[Continuous Data]\r\n');
warning('Setting uV as units, actual units unknown.');
%% data
for isamp = 1:size(indata.epochs.data,2)
    for ichan = 1:indata.num_sensors

        fprintf(fid,'%g',indata.epochs.data(ichan,isamp));
        if ichan<indata.num_sensors
            fprintf(fid,'\t');
        else
            fprintf(fid,'\r\n');
        end
    end
end
fclose(fid);
end
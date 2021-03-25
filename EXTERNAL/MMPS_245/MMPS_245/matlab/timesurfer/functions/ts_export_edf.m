function [HDR] = ts_export_edf(indata,outfname)
% [HDR] = ts_export_dat(indata,outfname)
% exports continuous data in timesurfer format to .edf file readable
% by neuroscan. 
% Required paramters:
% indata: Timesurfer datastruct containing continous data
% outfname: fullfile path to output file
%
% Created 02/20/13 BQR
% TODO: detect and accept average and epoch data inputs

%% Parse inputs
[Path,Name,Ext] = fileparts(outfname);
data = indata.epochs.data';

%% Confirm data Nsamps is divisable by sfreq
if mod(numel(indata.epochs.time),indata.sfreq)
   modsamps = floor(numel(indata.epochs.time)/indata.sfreq)*indata.sfreq;
   deltsamps = numel(indata.epochs.time)-modsamps;
   data = data(1:modsamps,:);
   warning('Data truncated by %i samples to ensure Nsamps is divisable Sfreq',deltsamps)
end


%% Assemble Header
HDR.FileName = outfname;
HDR.FILE.stdout = 1;
HDR.FILE.stderr = 2;
HDR.FILE.PERMISSION = 'w';
HDR.FILE.Path = Path;
HDR.FILE.Name = Name;
HDR.FILE.Ext = Ext(2:end);
HDR.FILE.size = [];
HDR.FILE.POS = 0;
HDR.TYPE = 'EDF';
% HDR.keycode
HDR.ERROR.status = 0;
HDR.ERROR.message = '';
HDR.NS = indata.num_sensors;
HDR.SampleRate = indata.sfreq;
HDR.T0 = clock;
HDR.Filter.Notch = nan(1,HDR.NS);
HDR.Filter.LowPass = nan(1,HDR.NS);
HDR.Filter.HighPass = nan(1,HDR.NS);
HDR.FLAG.FILT = 0;
HDR.FLAG.TRIGGERED = 0;
HDR.FLAG.UCAL = 0;
HDR.FLAG.OVERFLOWDETECTION = 1;
HDR.FLAG.FORCEALLCHANNEL = 0;
HDR.FLAG.OUTPUT = 'double';
HDR.EVENT.TYP = [];
HDR.EVENT.POS = [];
HDR.EVENT.CHN = [];
HDR.EVENT.DUR = [];
% HDR.ErroNo = [];
HDR.VERSION = 0;
HDR.Patient.Sex = 0;
HDR.Patient.Handedness = 'Unspecified';
HDR.Patient.Id = 0;
HDR.Patient.Name = '';
HDR.PID = 'Unspecified';
HDR.RID = pwd;
% HDR.Headlen = 25344;
HDR.reserved1 = 'Exported from Neuroscan SCAN 4 software     ';
% HDR.NRec = 437;
HDR.Dur = 1;
HDR.AS.H1 = '';
 HDR.AS.SPR =  repmat(HDR.SampleRate,HDR.NS,1);
%HDR.AS.SPR =  repmat(1,HDR.NS,1);
HDR.AS.SampleRate = repmat(HDR.SampleRate,HDR.NS,1);
HDR.AS.spb = HDR.SampleRate*HDR.NS;
HDR.AS.bi = (0:HDR.SampleRate:HDR.SampleRate*HDR.NS)';
HDR.AS.BPR = repmat(HDR.SampleRate*2,HDR.NS,1);
HDR.AS.SAMECHANTYP = 1;
HDR.AS.bpb = HDR.AS.spb*2;
HDR.AS.EVENTTABLEPOS = -1;
HDR.AS.c = HDR.SampleRate*HDR.NS;
HDR.AS.c2 = HDR.AS.bpb;
HDR.AS.TYP = 3;
HDR.Label = {indata.sensor_info.label};
HDR.Transducer = repmat({'Unknown. Usually EEG electrode.'},size(HDR.Label));
HDR.PhysDim = repmat({'uV'},size(HDR.Label));
warning('Setting uV as units, actual units unknown.');
HDR.PhysMax = repmat(5500,size(HDR.Label));
HDR.PhysMin = repmat(-5500,size(HDR.Label));
HDR.PhysMin = repmat(-5500,size(HDR.Label));
HDR.DigMin = repmat(-32768,size(HDR.Label));
HDR.DigMax = repmat(32767,size(HDR.Label));
HDR.PreFilt = repmat('HP:0.15 LP:30.00',HDR.NS,1);
HDR.GDFTYP = repmat(3,size(HDR.Label));
HDR.SPR = HDR.SampleRate;
HDR.THRESHOLD = [repmat(-32768,HDR.NS,1) repmat(32767,HDR.NS,1)];
% HDR.Cal = repmat(0.1678,size(HDR.Label));
% HDR.Off = repmat(0.0839,size(HDR.Label));
% HDR.Calib = 
HDR.PhysDimCode = repmat(4275,fliplr(size(HDR.Label)));
HDR.ELEC.XYZ = nan(HDR.NS,3);
HDR.ELEC.Phi = nan(HDR.NS,1);
HDR.ELEC.Theta = nan(HDR.NS,1);
HDR.LeadIdCode = nan(HDR.NS,1);
HDR.CHANTYP = repmat(' ',1,HDR.NS);
HDR.InChanSelect = (1:HDR.NS)';
HDR.SIE.RAW = 1;

%% manipulate paths
paths = path;
add_EEGLAB_flag = isempty(regexp(paths,'/usr/pubsw/packages/eeglab/eeglab9_0_5_6b','once'));
if add_EEGLAB_flag
    addpath(genpath('/usr/pubsw/packages/eeglab/eeglab9_0_5_6b'));
end

%% Export data
HDR = ts_ssave(HDR,data);

%% clean up paths
if add_EEGLAB_flag
    rmpath(genpath('/usr/pubsw/packages/eeglab/eeglab9_0_5_6b'));
end

function [HDR] = ts_ssave(FILENAME,DATA,TYPE,Fs,gdftyp)
% SSAVE saves signal data in various data formats
% 
% Currently are the following data formats supported: 
%    EDF, BDF, GDF, BKR, SND/AU, (WAV, AIF)
%    and WSCORE event file
%
% HDR = ssave(HDR,data);
% HDR = ssave(FILENAME,data,TYPE,Fs);
%
% FILENAME      name of file
% data  signal data, each column is a channel
% TYPE 	determines dataformat
% Fs	sampling rate	
%
% see also: SSAVE, SOPEN, SWRITE, SCLOSE, doc/README
%

% $Id: ssave.m 2205 2009-10-27 12:18:15Z schloegl $
% Copyright (C) 2003,2004,2007 by Alois Schloegl <a.schloegl@ieee.org>	
% This file is part of the biosig project http://biosig.sf.net/

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 2 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.



if isstruct(FILENAME),
        HDR = FILENAME;
        if isfield(HDR,'FileName'),
                FILENAME = HDR.FileName;
        else
                fprintf(2,'Error SSAVE: missing FileName.\n');	
                return; 
        end;
else
        HDR.FileName = FILENAME;
        HDR.SampleRate = Fs; 
        gdftyp = 16;
        HDR.TYPE = 'native';
end;

if (nargin > 1),
% 	[HDR.SPR, HDR.NS] = size(DATA); HDR.NRec = 1; 
%	HDR.AS = rmfield(HDR.AS,'SPR'); 
	if (strcmp(HDR.TYPE,'BDF') | strcmp(HDR.TYPE,'EDF') | strcmp(HDR.TYPE,'GDF')) & (~isfield(HDR,'DigMax') | ~isfield(HDR,'DigMin') |~isfield(HDR,'PhysMax') | ~isfield(HDR,'PhysMin'))
		HDR.PhysMax = max(DATA,[],1);
		HDR.PhysMin = min(DATA,[],1);
		ix = find(HDR.PhysMax == HDR.PhysMin);
		HDR.PhysMin(ix) = HDR.PhysMin(ix) - 1;
		if strcmp(HDR.TYPE,'BDF')
		   	[datatyp,HDR.THRESHOLD,datatypes,HDR.bits,HDR.GDFTYP] = gdfdatatype(511+24*ones(HDR.NS,1));
		   	HDR.DigMax = HDR.THRESHOLD(:,2)';
		   	HDR.DigMin = HDR.THRESHOLD(:,1)';
		elseif strcmp(HDR.TYPE,'EDF')
		   	[datatyp,HDR.THRESHOLD,datatypes,HDR.bits,HDR.GDFTYP] = gdfdatatype(3*ones(HDR.NS,1));
		   	HDR.DigMax = HDR.THRESHOLD(:,2)';
		   	HDR.DigMin = HDR.THRESHOLD(:,1)';
		elseif strcmp(HDR.TYPE,'GDF')
		   	[datatyp,HDR.THRESHOLD,datatypes,HDR.bits,HDR.GDFTYP] = gdfdatatype(16*ones(HDR.NS,1));
			HDR.DigMax = HDR.PhysMax;
			HDR.DigMin = HDR.PhysMin;
		else 	
			HDR.DigMax = HDR.PhysMax;
			HDR.DigMin = HDR.PhysMin;
		end; 
	end;    	

	if (nargin > 2),
        	if strcmp(TYPE,'GDF2'),
        		HDR.TYPE = 'GDF';
	        	HDR.VERSION = 2;
        	elseif strncmp(TYPE,'GDF',3),
        		HDR.TYPE = 'GDF';
      	 	 	HDR.VERSION = 1.25;
        	else	
	        	HDR.TYPE = TYPE; 	% type of data format
		end;        
	end;
	HDR = sopen(HDR,'w');
	HDR = swrite(HDR,DATA);
	HDR = sclose(HDR);
end;

% Convert EVENT into WSCORE event format
if all([length(HDR.EVENT.POS), length(HDR.EVENT.TYP)]),
	p = which('sopen'); [p,H,e] = fileparts(p);
	H = sload(fullfile(p,'../doc/eventcodes.txt'));

	HDR.EVENT.CodeDesc  = H.CodeDesc;
	HDR.EVENT.CodeIndex = H.CodeIndex;
	if isfield(HDR.EVENT,'DUR')
	        HDR.EVENT.POS = [HDR.EVENT.POS; HDR.EVENT.POS + HDR.EVENT.DUR];
	        HDR.EVENT.TYP = [HDR.EVENT.TYP; HDR.EVENT.TYP + hex2dec('8000')];
	end;
	OnOff = {'On','Off'};
	
	[HDR.EVENT.POS, ix] = sort(HDR.EVENT.POS);
	HDR.EVENT.TYP       = HDR.EVENT.TYP(ix);
	[TYP, IX, IY]       = unique(HDR.EVENT.TYP);

	% write "free form" scoring file for WSCORE
	fid   = fopen(fullfile(HDR.FILE.Path,[HDR.FILE.Name,'.C07']),'w');
	for k = 1:length(TYP), 
    		fprintf(fid,'%2i %s (%s)\r\n', k, HDR.EVENT.CodeDesc(mod(TYP(k),2^15)==HDR.EVENT.CodeIndex), OnOff{(TYP(k)>=2^15)+1});
	end;
	fclose(fid);

	% write "free form" scoring file for WSCORE
	fid = fopen(fullfile(HDR.FILE.Path,[HDR.FILE.Name,'.007']),'w');
	fprintf(fid,'%8i %i\r\n', [round(HDR.EVENT.POS(:)),IY(:)]');
	fclose(fid);
end;
end
end


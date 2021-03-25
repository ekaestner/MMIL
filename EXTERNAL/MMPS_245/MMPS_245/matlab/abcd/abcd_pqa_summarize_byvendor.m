function PQA_summarize_byvendor
%PQA_SUMMARIZE_BYSITE Summary of this function goes here
%   Detailed explanation goes here

[~,RootDirs] = MMIL_Get_ProjInfo('DAL_ABCD_PQA');

projroot = fileparts(fileparts(RootDirs.orig));
% indir - input directory containing jsonfile, also writes output
indir = sprintf('%s/ABCD_QA_proc/Summary',projroot);
outdir = strcat(indir,'/ByVendor');

if ~exist(outdir,'dir')
    cmd = sprintf('mkdir %s', outdir);
    unix(cmd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vendorList = {'GE';'SIEMENS';'Philips'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


listofcsv = dir(fullfile(indir,'*.csv'));

%%% Get last version

toRead = fullfile(indir,listofcsv(end).name);

summaryTable = readtable(toRead);

% ToDo: Function
%combinedScatterPlots(summaryTable, vendorList, outdir);

for i=1:length(vendorList)
    vendor = vendorList{i};
    v_Table = writeManufacturerTable(vendor,summaryTable,outdir);
    performancePlots(vendor,v_Table,outdir,1);
end

end

function [summaryTable] = writeManufacturerTable(vendor,summaryTable,outdir)

a = cellfun(@(x) strfind(x,vendor),summaryTable.Manufacturer, 'UniformOutput', false);
b = cellfun(@(x) isempty(x), a, 'UniformOutput', false);
toDelete = cell2mat(b);
summaryTable(toDelete,:) = [];

outdir = sprintf('%s/%s', outdir,vendor);
if ~exist(outdir,'dir')
    cmd = sprintf('mkdir %s', outdir);
    unix(cmd);
end

if ~isempty(summaryTable)
    outfile = sprintf('%s/%s.csv', outdir,vendor);
    writetable(summaryTable,outfile,'Delimiter',',','QuoteStrings',true);
end

end

function [multibandTable,nonmbTable_flip77,nonmbTable_flip10] = performancePlots(vendor, summaryTable, outdir, print)

%==========3 possible sequences============%

multiband_TR = 800;
conventional_TR = 2000;
flip_angle_77 = 77;
flip_angle_10 = 10;

multibandTable = summaryTable;
toKeep = summaryTable.RepetitionTime==multiband_TR;
multibandTable(~toKeep,:) = [];

nonmbTable = summaryTable;
toKeep = summaryTable.RepetitionTime==conventional_TR;
nonmbTable(~toKeep,:) = [];

nonmbTable_flip10 = nonmbTable;
toKeep = (nonmbTable.FlipAngle==flip_angle_10);
nonmbTable_flip10(~toKeep,:) = [];

nonmbTable_flip77 = nonmbTable;
toKeep = (nonmbTable.FlipAngle==flip_angle_77);
nonmbTable_flip77(~toKeep,:) = [];


if print
    if height(multibandTable)>1
        generatePlots(outdir, vendor, multibandTable, 'multiband-TR800');
    end

    if height(nonmbTable_flip77)>1
        generatePlots(outdir, vendor, nonmbTable_flip77, 'non-multiband-TR2000-flip77');
    end

    if height(nonmbTable_flip10)>1
        generatePlots(outdir, vendor, nonmbTable_flip10, 'non-multiband-TR2000-flip10');
    end
end


end


function generatePlots(outdir, vendor, subtable, name_append)

outdir = sprintf('%s/%s', outdir,vendor);
toPlot = {'Mean', 'StandardDev', 'SNR', 'SFNR', 'RMS', 'TemporalDrift', 'TemporalDrift_per_minute',...,
    'Max_Absolute_Temporal_Drift', 'RDC', 'FWHM_x', 'FWHM_y', 'FWHM_z', 'MeanGhost', 'TopGhost',...,
    'ImageFrequency'};

[subtable] = sortrows(subtable,{'StudyDate'},{'ascend'});

for i=1:height(subtable)
    fulldate = cell2mat(subtable.StudyDate(i));
    year = str2double(fulldate(1:4));
    month = str2double(fulldate(5:6));
    day = str2double(fulldate(7:8));
    ndate(i) = datenum(year,month,day);
end

for i=1:length(toPlot)
    index_table = find(strcmpi(subtable.Properties.VariableNames,toPlot(i)));
    y = subtable{:,index_table};
    plot(ndate,y, '-o');grid;
    datetick('x',22);
    set(gca,'XMinorTick','on');
    limits=[0.5*mean(y) 1.5*mean(y), min(y), max(y)];
    ylim([min(limits) max(limits)]);
    
    title(sprintf('%s-%s',vendor,name_append));
    xlabel('Date');
    ylabel(sprintf('%s',cell2mat(toPlot(i))));

    set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')'); % this disables the figure and set the 'CreateFcn' property simultaneously
    resultFile = sprintf('%s_%s_%s.tiff',vendor,cell2mat(toPlot(i)),name_append);
    file = fullfile(outdir, resultFile);
    saveas(gcf, file) % save the figure as a .tiff folder
    close;    

end
end


%TO BE COMPLETED (call commented in the code)

function combinedScatterPlots(summaryTable, vendorList, outdir)

ge_Table = summaryTable;
a = cellfun(@(x) strfind(x,vendorList{1}),summaryTable.Manufacturer, 'UniformOutput', false);
b = cellfun(@(x) isempty(x), a, 'UniformOutput', false);
toDelete = cell2mat(b);
ge_Table(toDelete,:) = [];

siemens_Table = summaryTable;
a = cellfun(@(x) strfind(x,vendorList{2}),summaryTable.Manufacturer, 'UniformOutput', false);
b = cellfun(@(x) isempty(x), a, 'UniformOutput', false);
toDelete = cell2mat(b);
siemens_Table(toDelete,:) = [];

philips_Table = summaryTable;
a = cellfun(@(x) strfind(x,vendorList{3}),summaryTable.Manufacturer, 'UniformOutput', false);
b = cellfun(@(x) isempty(x), a, 'UniformOutput', false);
toDelete = cell2mat(b);
philips_Table(toDelete,:) = [];

[multibandTable_ge,nonmbTable_flip77_ge,nonmbTable_flip10_ge] = performancePlots('', ge_Table, '', 0);
[multibandTable_siemens,nonmbTable_flip77_siemens,nonmbTable_flip10_siemens] = performancePlots('', siemens_Table, '', 0);
[multibandTable_philips,nonmbTable_flip77_philips,nonmbTable_flip10_philips] = performancePlots('', philips_Table, '', 0);


toPlot = {'Mean', 'StandardDev', 'SNR', 'SFNR', 'RMS', 'TemporalDrift', 'TemporalDrift_per_minute',...,
    'Max_Absolute_Temporal_Drift', 'RDC', 'FWHM_x', 'FWHM_y', 'FWHM_z', 'MeanGhost', 'TopGhost',...,
    'ImageFrequency'};

for i=1:length(toPlot)
    index_table = find(strcmpi(multibandTable_ge.Properties.VariableNames,toPlot(i)));
    ge_param = multibandTable_ge{:,index_table};
    index_table = find(strcmpi(multibandTable_siemens.Properties.VariableNames,toPlot(i)));
    siemens_param = multibandTable_ge{:,index_table};
    index_table = find(strcmpi(multibandTable_philips.Properties.VariableNames,toPlot(i)));
    philips_param = multibandTable_ge{:,index_table};
    
    compound_vals = [ge_param; siemens_param; philips_param];
    compound_indexes = [ones(length(ge_param),1);2*ones(length(siemens_param),1);3*ones(length(ge_param),1)];
    
    gscatter(compound_indexes, compound_vals, compound_indexes);
end


end



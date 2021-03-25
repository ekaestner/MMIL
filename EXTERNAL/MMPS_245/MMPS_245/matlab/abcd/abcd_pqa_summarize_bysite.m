function PQA_summarize_bysite
%PQA_SUMMARIZE_BYSITE Summary of this function goes here
%   Detailed explanation goes here
[~,RootDirs] = MMIL_Get_ProjInfo('DAL_ABCD_PQA');

projroot = fileparts(fileparts(RootDirs.orig));
% indir - input directory containing jsonfile, also writes output
indir = sprintf('%s/ABCD_QA_proc/Summary',projroot);
outdir = strcat(indir,'/BySite');

if ~exist(outdir,'dir')
    cmd = sprintf('mkdir %s', outdir);
    unix(cmd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%siteList = {'SRI'; 'CHLA';'UMN';'UCSD';'LIBR';'UMICH';'UPMC';'UTAH';'CUB';...,
%    'WUSTL';'FIU';'MSSM';'OHSU';'UCLA';'UFL';'OAHU';'UVM';'VCU';'MUSC';'UWM';'YALE'};

siteList = {'SRI';'CHLA';'UMN_1';'UMN_2';'UCSD_1';'UCSD_2';'LIBR_1';'LIBR_2';'UMICH_1';'UMICH_2';'UPMC';'UTAH_1';'UTAH_2';'CUB';...
            'WUSTL_1';'WUSTL_2';'FIU';'MSSM';'OHSU';'UCLA';'UFL';'UMB_1';'UMB_2';'UVM';'VCU';'MUSC';'UWM';'YALE_1';'YALE_2';'ROC'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


listofcsv = dir(fullfile(indir,'*.csv'));

%%% Get last version

toRead = fullfile(indir,listofcsv(end).name);

summaryTable = readtable(toRead);

for i=1:length(siteList)
    site = siteList{i};
    writeSiteTable(site,summaryTable,outdir);
    performancePlots(site,summaryTable,outdir);
end

end

function writeSiteTable(site,summaryTable,outdir)

toDelete = ~(strcmp(summaryTable.Site, site));
summaryTable(toDelete,:) = [];

outdir = sprintf('%s/%s', outdir,site);
if ~exist(outdir,'dir')
    cmd = sprintf('mkdir %s', outdir);
    unix(cmd);
end

if ~isempty(summaryTable)
    outfile = sprintf('%s/%s.csv', outdir,site);
    writetable(summaryTable,outfile,'Delimiter',',','QuoteStrings',true);
end

end

function performancePlots(site, summaryTable, outdir)

%==========3 possible sequences============%

toDelete = ~(strcmp(summaryTable.Site, site));
summaryTable(toDelete,:) = [];

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

if height(multibandTable)>1
    generatePlots(outdir, site, multibandTable, 'multiband_TR800');
end

if height(nonmbTable_flip77)>1
    generatePlots(outdir, site, nonmbTable_flip77, 'nonmultiband_TR2000_flip77');
end

if height(nonmbTable_flip10)>1
    generatePlots(outdir, site, nonmbTable_flip10, 'nonmultiband_TR2000_flip10');
end



end


function generatePlots(outdir, site, subtable, name_append)

outdir = sprintf('%s/%s', outdir,site);
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
    
    title(sprintf('%s-%s',site,name_append));
    xlabel('Date');
    ylabel(sprintf('%s',cell2mat(toPlot(i))));

    set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')'); % this disables the figure and set the 'CreateFcn' property simultaneously
    resultFile = sprintf('%s_%s_%s.tiff',site,cell2mat(toPlot(i)),name_append);
    file = fullfile(outdir, resultFile);
    saveas(gcf, file) % save the figure as a .tiff folder
    close;
        

end

end


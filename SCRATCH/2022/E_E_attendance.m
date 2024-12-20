clear; clc;

atn_fld = '/home/ekaestner/Dropbox/McDonald Lab/Donatello/enigma/calls/Attendances';

%% Load attendances
atn_fle = dir(atn_fld);
    atn_fle = {atn_fle(3:end-1).name};
    
atn_dta     = cell(1,numel(atn_fle));
atn_dta_tme = nan(1,numel(atn_fle));
atn_dta_col = cellfun(@(x) x(1:end-4),atn_fle,'uni',0);
mss_ind     = zeros(1,numel(atn_fle));
atn_dta_typ = cellfun(@(x) x(end),regexp(atn_dta_col,'_','split'));
atn_dta_yer = cellfun(@str2num,cellfun(@(x) x(1),regexp(atn_dta_col,'_','split')));

for iF = 1:numel(atn_fle)
    try 
        atn_dta{iF} = mmil_readtext([ atn_fld '/' atn_fle{iF} ]);
        atn_dta_tme(iF) = atn_dta{iF}{1};
        atn_dta{iF}     = atn_dta{iF}(2:end);
    catch
        atn_dta{iF} = '';
        atn_dta_tme(iF) = NaN;
        mss_ind(iF) = 1;
    end
end

cell2csv( [ atn_fld '/' 'characterization' '/' '01_missing_meetings.txt' ], atn_dta_col(logical(mss_ind))' );

%% Get total number of attendees
atn_dta_num = cellfun(@numel,atn_dta);
    atn_dta_num( atn_dta_num==0 ) = NaN;

% Plot by time
fcfg = [];

fcfg.xdt     = { atn_dta_tme(strcmpi(atn_dta_typ,'TC'))' atn_dta_tme(strcmpi(atn_dta_typ,'Methods'))' atn_dta_tme(strcmpi(atn_dta_typ,'Cognitive'))'  };
fcfg.ydt     = { atn_dta_num(strcmpi(atn_dta_typ,'TC'))' atn_dta_num(strcmpi(atn_dta_typ,'Methods'))' atn_dta_num(strcmpi(atn_dta_typ,'Cognitive'))'  };

fcfg.fce_col = {rgb('green') rgb('purple') rgb('burnt orange') };
fcfg.edg_col = repmat({rgb('black')},1,numel(fcfg.xdt));

fcfg.xlb = { 'Time of Day (PST)' };
fcfg.ylb = { '# Attendees'  };

fcfg.jtr = 1;
%fcfg.xlm = [ -5 5 ];
fcfg.ylm = [ 0 max(cat(1,fcfg.ydt{:}))+5 ];

fcfg.out_dir = [ atn_fld '/' 'characterization' '/' ];
fcfg.out_nme = '02_Number_of_attendees_by_time';

ejk_scatter(fcfg)

% Plot by year
fcfg = [];

fcfg.xdt     = { atn_dta_yer(strcmpi(atn_dta_typ,'TC'))' atn_dta_yer(strcmpi(atn_dta_typ,'Methods'))' atn_dta_yer(strcmpi(atn_dta_typ,'Cognitive'))'  };
fcfg.ydt     = { atn_dta_num(strcmpi(atn_dta_typ,'TC'))' atn_dta_num(strcmpi(atn_dta_typ,'Methods'))' atn_dta_num(strcmpi(atn_dta_typ,'Cognitive'))'  };

fcfg.fce_col = {rgb('green') rgb('purple') rgb('burnt orange') };
fcfg.edg_col = repmat({rgb('black')},1,numel(fcfg.xdt));

fcfg.xlb = { 'Time of Day (PST)' };
fcfg.ylb = { '# Attendees'  };

fcfg.jtr = 1;
fcfg.xlm = [ 2015 2023 ];
fcfg.ylm = [ 0 max(cat(1,fcfg.ydt{:}))+5 ];

fcfg.out_dir = [ atn_fld '/' 'characterization' '/' ];
fcfg.out_nme = '03_Number_of_attendees_by_year';

ejk_scatter(fcfg)

%% Put together list of attendees 
tot_nme = cat(1,atn_dta{:});

% Who is attending
tot_nme_tbl = tabulate(tot_nme);
    [~, tot_nme_srt] = sort(cell2mat(tot_nme_tbl(:,2)));
    tot_nme_tbl = tot_nme_tbl(flipud(tot_nme_srt),1:2);

cell2csv( [ atn_fld '/' 'characterization' '/' '04_top_attendees.txt' ], tot_nme_tbl )
    
% When are they attending
tot_atn_tbl = [tot_nme_tbl(:,1) repmat({''},size(tot_nme_tbl,1),numel(atn_fle))];
for iC = 1:numel(atn_fle)
    for iS = 1:size(tot_atn_tbl,1)
        if any(strcmp(atn_dta{iC},tot_atn_tbl{iS,1}))
            tot_atn_tbl{iS,iC+1} = 'Y';
        end
    end
end

cell2csv( [ atn_fld '/' 'characterization' '/' '05_attendance_record.txt' ], [ {'name'} atn_fle ;  tot_atn_tbl ] );


function [pGUID,EventName,SessionType] = abcd_get_pGUID(fstem)
%function [pGUID,EventName,SessionType] = abcd_get_pGUID(fstem)
%
% Purpose: extract pGUID and EventName from fstem
%
% Input:
%   fstem: string from json/tgz file stem
%
% Output:
%   pGUID: globally unique ID
%   EventName: name of data acquisition "event"
%     e.g. 'baseline_year_1_arm_1'
%   SessionType: type of data acquisition session
%     e.g. 'SessionA1'
%
% Created:  09/13/16 by Don Hagler
% Last Mod: 09/16/16 by Don Hagler
%

pGUID = []; EventName = []; SessionType = [];

if ~mmil_check_nargs(nargin,1), return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

event_list = {...
  'screener_arm_1'...
  'baseline_year_1_arm_1'...
  '6_month_follow_up_arm_1'...
  '1_year_follow_up_y_arm_1'...
  '18_month_follow_up_arm_1'...
  '2_year_follow_up_y_arm_1'...
  '30_month_follow_up_arm_1'...
  '3_year_follow_up_y_arm_1'...
  '42_month_follow_up_arm_1'...
  '4_year_follow_up_y_arm_1'...
  '54_month_follow_up_arm_1'...
  '5_year_follow_up_y_arm_1'...
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(event_list)
  event_name = event_list{i};
  pat = sprintf('(?<ID>.+)_%s_(?<Sess>Session[^_]+)_',event_name);
  k = regexp(fstem,pat,'names');
  if ~isempty(k)
    pGUID = k.ID;
    EventName = event_name;
    SessionType =  k.Sess;
    break;
  end;
end;


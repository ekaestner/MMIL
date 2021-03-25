function gain_chan = ts_gain_intpnt2chan(coil_info,gain_intpnt)
%function gain_chan = ts_gain_intpnt2chan(coil_info,gain_intpnt)
%
% Purpose: transform gain matrix from integration points to channels
%
% Required Input:
%   coil_info: struct array with n_chan members
%     coil_info(i).n: number of integration point for chan i
%     coil_info(i).wei: diag weighting matrix for chan i
%   gain_intpnt: gain matrix with size = [n_intpnt,n_sources]
%
% Output:
%   gain_chan: gain matrix with size = [n_chan,n_sources]
%
% Note: based on gain_intpnt2chan by M.X. Huang
%
% Created:  10/11/13 by Don Hagler
% Last Mod: 10/11/13 by Don Hagler
%

n_chan=length(coil_info);
index=0;
n_col=size(gain_intpnt,2);
gain_chan=zeros(n_chan,n_col);

for i=1:n_chan
  n_intpnt=coil_info(i).n;
  wei=coil_info(i).wei;
  index_intpnt=index+[1:n_intpnt];
  gain_chan(i,:)=sum(wei*gain_intpnt(index_intpnt,:),1);
  index=index+n_intpnt;
end    


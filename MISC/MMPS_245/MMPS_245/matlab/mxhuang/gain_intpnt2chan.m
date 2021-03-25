function gain_chan=gain_intpnt2chan(coil_info,gain_intpnt);
% function gain_chan=gain_intpnt2chan(coil_info,gain_intpnt);
% transform gain matrix from integration points to channels
%
% coil_info: 1xn_chan structure. 
%   coil_info(i).n: number of integration point for chann i
%   coil_info(i).wei: diag weighting matrix for chann i

n_chann=length(coil_info);
index=0;
n_col=size(gain_intpnt,2);
gain_chan=zeros(n_chann,n_col);

for i=1:n_chann
    n_intpnt=coil_info(i).n;
    wei=coil_info(i).wei;
    index_intpnt=index+[1:n_intpnt];
    gain_chan(i,:)=sum(wei*gain_intpnt(index_intpnt,:),1);
    index=index+n_intpnt;
end    
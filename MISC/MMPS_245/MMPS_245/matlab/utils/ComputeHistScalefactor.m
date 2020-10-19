function sf = ComputeHistScalefactor(hc,bv,hc_template,bv_template,bn0_template)

v1 = cumsum(hc_template(bn0_template:end))/sum(hc_template(bn0_template:end));
flist = exp([-1:0.01:1]);
costvec = zeros(1,length(flist));
for fi = 1:length(flist)
  f = flist(fi);
  ind = [1:length(hc_template)]*f;
  ind = max(1,min(length(hc_template),ind));
  hc_res = interp1(hc,ind);
  v2_res = cumsum(hc_res(bn0_template:end))/sum(hc_res(bn0_template:end));
  costvec(fi) = sqrt(mean((v1-v2_res).^2));
%  figure(10); plot([v1' v2_res']); fprintf(1,'fi=%d, f=%f\n',fi,f); pause;
end

%figure; plot(flist,costvec)

[minv,mini] = min(costvec);
sf = bv_template(end)/(flist(mini)*bv(end));

%keyboard

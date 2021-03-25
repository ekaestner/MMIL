function ctx_vol = ctx_load_cor(fname)

[vol,M] = load_cor(fname);
ctx_vol = ctx_mgh2ctx(vol,M);


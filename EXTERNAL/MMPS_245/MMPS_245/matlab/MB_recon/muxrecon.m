function muxrecon(pfile, outfile, ext_cal, slice, n_vcoils, recon_method)

fprintf('mux_epi_main("%s", "%s", "%s", %d, [], %d, 0, "%s", 1, 1, 0)\n', pfile, outfile, ext_cal, str2double(slice), str2double(n_vcoils), recon_method);
mux_epi_main(pfile, outfile, ext_cal, str2double(slice), [], str2double(n_vcoils), 0, recon_method, 1, 1, 0);

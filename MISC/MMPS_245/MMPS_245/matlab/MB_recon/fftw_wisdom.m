

% Estimate fftw wisdom for a bunch of common mux recon sizes. For some
% matrix sizes, the speed-up over the conventional heuristics can be 25% or
% more.
fftw('planner','exhaustive');

sz = [64,170]
dat = rand(sz(2), sz(2));
for ii=sz(1):sz(2)
    fprintf('Computing %d of %d ffts...\n', ii-sz(1)+1, sz(2)-sz(1));
    fflush(stdout);
    ft = ifft2(dat(1:ii,1:ii));
end

wisdom = fftw('dwisdom');
save('fftw_wisdom.txt', 'wisdom');

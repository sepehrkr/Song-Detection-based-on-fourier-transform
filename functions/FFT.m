function single_sided_power_spectrum = FFT(X)
    %%% single_sided_power_spectrum
    L = length(X);
    Y = fft(X);
    two_side_spec = (abs(Y)/L) .^ 2;
    End = floor(L/2)+1;
    single_sided_power_spectrum = two_side_spec(1:End) .* 2;
    single_sided_power_spectrum(1) = single_sided_power_spectrum(1)/2; 
end
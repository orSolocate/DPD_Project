function new_sig = resampling(sig_in, old_freq, new_freq)

    [P, Q] = rat(new_freq/old_freq);
    new_sig = resample(sig_in,P,Q);

end
clear all
clf;
load DSPI_StereoRadio_MATLAB.mat
Fs = 400e3;
Ts = 1/(Fs);
Fs_mic = 16e3;

RXw =  fft(RXn);
num_samples = length(RXw);
freq_bin_factor = Fs/num_samples;
freq = freq_bin_factor*(0:num_samples-1)';

mono_raw = modulate_signal(RXw,70e3,Fs);
mono_fft = ideal_lowpass(mono_raw, 5e3,Fs);
mono = ifft(mono_fft);
mono_hat = downsample(mono, Fs/Fs_mic);


sub_raw = modulate_signal(RXw,90e3,Fs);
sub_fft = ideal_lowpass(sub_raw,5e3,Fs);
sub = ifft(sub_fft);
sub_hat = downsample(sub, Fs/Fs_mic);

num_mic_samples = length(mono_hat);
freq_mic = Fs_mic/num_mic_samples*(0:num_mic_samples-1)';

left_hat = 0.5*(mono_hat+sub_hat);
right_hat = 0.5*(mono_hat-sub_hat);

sound(2*right_hat,Fs_mic)
input("Press enter")
sound(2*left_hat,Fs_mic)

function filtered_sig = ideal_lowpass(signal_fft,cutoff_freq,Fs)
    num_samples = length(signal_fft);
    passband_freq_index = find_freq_index(cutoff_freq,Fs, num_samples);
    rectangle = zeros(size(signal_fft));
    rectangle(1:passband_freq_index+1) = 1;
    rectangle(end-passband_freq_index+1:end) = 1;
    filtered_sig = rectangle .* signal_fft;
end

function modsig = modulate_signal(signal_fft, shift_down_freq, Fs) 
    time_sig = ifft(signal_fft);
    num_samples = length(signal_fft);
    mod_time_sig = time_sig .* cos(2*pi*shift_down_freq*(0:num_samples-1)/Fs)';
    modsig = fft(mod_time_sig);
end

function index = find_freq_index(freq,Fs, N)
    index = floor(freq*N/Fs);
end
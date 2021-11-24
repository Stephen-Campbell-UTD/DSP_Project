clear all
clf;
load DSPI_StereoRadio_MATLAB.mat
Fs = 400e3;
Ts = 1/(Fs);


RXw =  fft(RXn);
num_samples = length(RXw);
freq_bin_factor = Fs/num_samples;
freq = freq_bin_factor*(0:num_samples-1)';
mono_raw = modulate_signal(RXw,70e3,Fs);
mono_unfiltered = ifft(mono_raw);
mono = lowpass(mono_unfiltered, 5e3, Fs, 'ImpulseResponse','fir', 'Steepness', 0.99);
mono_fft = fft(mono);


sub_raw = modulate_signal(RXw,90e3,Fs);
sub_unfiltered = ifft(sub_raw);
sub = lowpass(sub_unfiltered, 5e3, Fs, 'ImpulseResponse','fir', 'Steepness', 0.99);
sub_fft = fft(sub);

subplot(3,1,1)
plot(freq,abs(RXw))
subplot(3,1,2)
hold on
plot(freq,abs(mono_raw), 'r')
plot(freq,abs(mono_fft), 'b')
subplot(3,1,3)
hold on
plot(freq,abs(sub_raw) ,'r')
plot(freq,abs(sub_fft), 'b')
hold off

dec_mono = downsample(mono,2);
dec_sub = downsample(sub,2);

% sound(dec_sub, Fs/2)
% figure
% plot(dec_mono)
% sound(2*dec_mono,Fs/2)
left = 0.5*(mono+sub);
right = 0.5*(mono-sub);
dec_left = downsample(left,2);
dec_right = downsample(right,2);
figure
subplot(2,1,1);
plot(left)
subplot(2,1,2);
plot(right)
sound(5*dec_left,Fs/2)


function filtered_sig = ideal_lowpass(signal_fft,cutoff_freq,Fs)
    num_samples = length(signal_fft);
    passband_freq_index = find_freq_index(cutoff_freq,Fs, num_samples);
    rectangle = zeros(size(signal_fft));
    rectangle(1:passband_freq_index+1) = 1;
    rectangle(end-passband_freq_index+1:end) = 1;
    filtered_sig = rectangle .* signal_fft;
end

function modsig = modulate_signal(signal_fft, shift_down_freq, Fs)
    
    num_samples = length(signal_fft);
    shift_freq_index = find_freq_index(shift_down_freq,Fs,num_samples);
    first_term = circshift(signal_fft,-shift_freq_index);
    second_term = circshift(signal_fft,-(mod(-shift_freq_index,num_samples)));
    modsig = first_term + second_term;
    
%     time_sig = ifft(signal_fft);
%     num_samples = length(signal_fft);
%     mod_time_sig = time_sig .* cos(2*pi*shift_down_freq*(1:num_samples)/Fs)';
%     modsig = fft(mod_time_sig);
end

function index = find_freq_index(freq,Fs, N)
    index = floor(freq*N/Fs);
%     index = freq*
end

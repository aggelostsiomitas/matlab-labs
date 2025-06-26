
clear;
clc;
signal = load('signal.txt');
signal = signal(:)';
fs = 1000; 
N = length(signal); 
t = (0:N-1)/fs; 


frequencies = 1:0.1:350; 
cosine_basis = cos(2 * pi * frequencies' * t);
inner_products = abs(cosine_basis * signal');

[~, peak_idx] = maxk(inner_products, 4);
estimated_frequencies = frequencies(peak_idx);


reconstructed_signal = zeros(size(t));
for f = estimated_frequencies
    reconstructed_signal = reconstructed_signal + cos(2 * pi * f * t);
end


fine_frequencies = [];
for f = estimated_frequencies
    fine_frequencies = [fine_frequencies, f-0.5:0.0001:f+0.5];
end
frequencies2 = unique(fine_frequencies);

new_cosine_basis = cos(2 * pi * frequencies2' * t);
inner_products2 = abs(new_cosine_basis * signal');

[~, peak_idx2] = maxk(inner_products2, 4);
new_estimated_frequencies = frequencies2(peak_idx2);

new_reconstructed_signal = zeros(size(t));
for f2 = new_estimated_frequencies
    new_reconstructed_signal = new_reconstructed_signal + cos(2 * pi * f2 * t);
end


figure;
subplot(2,2,1);
plot(t, signal, 'r', t, reconstructed_signal, 'g--');
title('Original vs. Reconstructed Signal (Coarse Resolution)');
legend('Original Signal', 'Reconstructed Signal');

subplot(2,2,2);
plot(t, abs(signal - reconstructed_signal));
title('Absolute Difference (Coarse Resolution)');

subplot(2,2,3);
plot(t, signal, 'r', t, new_reconstructed_signal, 'b--');
title('Original vs. Reconstructed Signal (High Resolution)');
legend('Original Signal', 'Reconstructed Signal (High Resolution)');

subplot(2,2,4);
plot(t, abs(signal - new_reconstructed_signal));
title('Absolute Difference (High Resolution)');

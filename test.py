from iterative_fft import fft, nearest_pow_2
import librosa
import numpy as np

original_audio, sample_rate = librosa.load('audio_files/River Flows in You.wav', sr=None)

reconstructed_audio = original_audio

# Applying FFT and IFFT consecutively to test for accuracy issues with compounding errors
number_of_iters = 20
for i in range(number_of_iters):
    reconstructed_audio = fft(fft(reconstructed_audio, 0), 1)

assert(np.allclose(np.append(original_audio, [0]*(nearest_pow_2(len(original_audio))-len(original_audio))), np.array(reconstructed_audio), atol=1e-6))
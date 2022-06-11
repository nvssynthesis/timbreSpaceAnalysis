# timbreSpaceAnalysis

load an audio file or folder of audio files to be analyzed for pitch, loudness, and timbral features.
outputs a .yaml file.
the intended use is to be used to train a neural network (i provide one at https://github.com/nvssynthesis/essentia2keras) to map these features to spectra.
the spectra contained in the created .yaml file are formatted specially to be readily available for wavetable generation.
the intended wavetable oscillator that uses this format is wtianns~ (https://github.com/nvssynthesis/wtianns).

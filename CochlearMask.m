% function [cochlearWeights, centreFreq, cochlearBandIdx, cochlearBinWidth] = ...
%          CochlearMask(nFreqBins, sFreq, nChan, loFreq, hiFreq, spl) 
%
% CochlearMask calculates a cochlear STFT Mask for cochlear-like analyses
%
%%%%%%%%%%%%%%%%%%%% Output Arguments %%%%%%%%%%%%%%%%%%%%
% cochlearWeights are the FFT Weights for each cochlear channel
% cochlearBandIdx are the FFT bins corresponding to the end of each
%                 cochlear channel
% cochlearBinWidth is the Bin Width for each cochlear channel
%%%%%%%%%%%%%%%%%%%%%%% Mandatory Arguments %%%%%%%%%%%%%%%%
% nFreqBins is the number of FFT bins to use
% sFreq is the sample frequency
% nChan is the number of cochlear channels
% loFreq is the low frequency end for the cochlear channels
% hiFreq is the high frequency end for the cochlear channels
% spl is the sound level in dB SPL

function [cochlearWeight, centreFreq, cochlearBandIdx, cochlearBinWidth] = CochlearMask(nFreqBins, sFreq, nChan, loFreq, hiFreq, spl) 
      
%%%% Fill in the function


% function [audfilt,audfreq] = roex(freqhz,spl,lowfreq,highfreq,Nfreqbins,sfreq);
%
% roex calculates a roex auditory filter centered at freqhz.
% The algorithm is the one specified in:
% Derivation of auditory filter shapes from notched-noise data
% by Brian Glasberg and Brian Moore, Hearing Research vol. 47,
% pp. 103-138, 1990 and also:
% Formulae describing frequency selectivity as a function of
% frequency and level, and their use in calculating excitation patterns
% by Brian Moore and Brian Glasberg, Hearing Research vol. 28,
% pp. 209-225, 1987. 
% Note that the low-frequency skirt of the filter changes shape
% with dB SPL.
%
%%%%%%%%%%%%%%%%%%%% Output Arguments %%%%%%%%%%%%%%%%%%%%
% audfilt is a roex cochlear filter at freqhz
% audfreq is the frequency axis for this filter
%%%%%%%%%%%%%%%%%%%%%%% Mandatory Arguments %%%%%%%%%%%%%%%%
% freqhz is the frequency in Hz
% spl is the sound level in dB SPL
% lowfreq is the low frequency limit at which to evaluate the filter 
% highfreq is the high frequency limit at which to evaluate the filter 
% Nfreqbins is the number of FFT points
function [audfilt,audfreq] = roex(freqhz,spl,lowfreq,highfreq,Nfreqbins,sfreq);

% Calculate frequency indices
hfindex = floor((highfreq*Nfreqbins/sfreq)+1);
lfindex = ceil((lowfreq*Nfreqbins/sfreq)+1);

% Calculate index at freqhz
freqindex = round((freqhz*Nfreqbins/sfreq)+1);

% Make full freq vector
freq = ([0:Nfreqbins-1]/Nfreqbins)*sfreq;

% Setup Constants
c(1) = 24.673;
c(2) = 4.368;
c(3) = 2302.6 / (c(1) * c(2));

% Calc Bandwidth at freqhz
ERB = c(1)*(c(2)*(freqhz/1000) + 1.0);

% Calc dB SPL/ERB
dBsplerb = spl;

% Calc p51 numbers
p51_1k = 4.0*1000.0/(c(1)*(c(2)+1.0));
p51 = 4.0*freqhz/ERB;

low_g = abs((freq(lfindex:freqindex) - freqhz)/freqhz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%p_lo = p51 - 0.380*(p51/p51_1k)*(dBsplerb - 51.0);
% A change in their formula, see
% A model for the prediction of thresholds, loudness
% and partial loudness, J. Audio Eng. Soc. Vol 45, No 4, 1997,
% page 224-240.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_lo = p51 - 0.35*(p51/p51_1k)*(dBsplerb - 51.0);
lofilt = (1.0 + p_lo*low_g).*exp(-p_lo*low_g);

hi_g = abs((freq(freqindex+1:hfindex) - freqhz)/freqhz);
p_hi = p51;
hifilt = (1.0 + p_hi*hi_g).*exp(-p_hi*hi_g);
audfilt = [lofilt,hifilt];
audfreq = freq(lfindex:hfindex);


cochlearWeight = audfilt;
centreFreq = audfreq;
cochlearBandIdx = freq;
cochlearBinWidth = ERB;


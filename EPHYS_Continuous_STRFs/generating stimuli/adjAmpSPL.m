function [Stim, d] = adjAmpSPL(Stim,Amp,calib_level,calib_V)

sigRMS = rms(Stim);
sigVolts = sigRMS*sqrt(2);
a = 10^((Amp-calib_level)/20)*calib_V;
Stim = Stim*a/sigVolts;
d = 20*log10(a/sigVolts); % dB value added to stim to get to Amp level
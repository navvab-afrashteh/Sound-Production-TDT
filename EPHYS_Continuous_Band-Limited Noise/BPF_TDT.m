function Hd = BPF_TDT(Fs, N, Fpass1, Fpass2)
%BPF_TDT Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 8.2 and the Signal Processing Toolbox 6.20.
% Generated on: 21-Mar-2018 16:33:55

% FIR least-squares Bandpass filter designed using the FIRLS function.

% All frequency values are in Hz.
Fstop1 = Fpass1*0.95;     % First Stopband Frequency
Fstop2 = Fpass2*1.02;    % Second Stopband Frequency
Wstop1 = 1;       % First Stopband Weight
Wpass  = 0.1;     % Passband Weight
Wstop2 = 1;       % Second Stopband Weight

% Calculate the coefficients using the FIRLS function.
b  = firls(N, [0 Fstop1 Fpass1 Fpass2 Fstop2 Fs/2]/(Fs/2), [0 0 1 1 0 ...
           0], [Wstop1 Wpass Wstop2]);
Hd = dfilt.dffir(b);

% [EOF]

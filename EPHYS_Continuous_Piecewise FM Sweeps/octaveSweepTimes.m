function SweepTimes = octaveSweepTimes(Pulse_Dur_spec, Log_Levels)

T1 = Pulse_Dur_spec(1); % the minimum pulse duration (linear scale, msec)
T2 = Pulse_Dur_spec(end-1); % the maxmimum pulse duration  (linear scale, kHz)
Tpoints = Pulse_Dur_spec(end); % number of pulse durations
if Log_Levels
    SweepTimes = 2.^linspace(log2(T1),log2(T2),Tpoints);
else
    SweepTimes = linspace(T1,T2,Tpoints);
end

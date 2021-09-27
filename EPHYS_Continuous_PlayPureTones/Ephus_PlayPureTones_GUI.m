clc
clear
close all

puretonesGUI = figure('units','normalized','position',[0.2 0.5 0.45 0.45],...
    'menubar','none','numbertitle','off','name','Pure Tones Specifications');

xtext = 0.00;
Ltext = 0.39;
xval = xtext + Ltext;
Lval = 0.15;

min_Freq_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','Minimum Frequency (kHz) = ', 'units','normalized',...
    'position',[xtext 0.9 Ltext 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'horizontalalignment','right');

min_Freq_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string',4,'units','normalized','position',[xval 0.9 Lval 0.065],...
    'tag','min_Freq');


max_Freq_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','Maximum Frequency (kHz) = ', 'units','normalized',...
    'position',[xtext 0.8 Ltext 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'horizontalalignment','right');

max_Freq_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string',32,'units','normalized','position',[xval 0.8 Lval 0.065],...
    'tag','max_Freq');


Freq_per_Oct_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','Frequency per Octave = ', 'units','normalized',...
    'position',[xtext 0.7 Ltext 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'horizontalalignment','right');

Freq_per_Oct_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string',0.667,'units','normalized','position',[xval 0.7 Lval 0.065],...
    'tag','Freq_per_Oct');


Amp_Array_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','Amplitude Array (dB) = ', 'units','normalized',...
    'position',[xtext 0.6 Ltext 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'horizontalalignment','right');

Amp_Array_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string','[40:20:80]','units','normalized','position',[xval 0.6 1.0*Lval 0.065],...
    'fontsize',7,'tag','Amp_Array');


ramp_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string',['ramp (msec)(change TDT) = '], 'units','normalized',...
    'position',[xtext 0.5 Ltext 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
	'Max',2,'Min',0,...
    'horizontalalignment','right');

ramp_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string',5,'units','normalized','position',[xval 0.5 Lval 0.065],...
    'tag','ramp');


Pulse_Dur_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','Pulse Duration (msec) = ', 'units','normalized',...
    'position',[xtext 0.4 Ltext 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'horizontalalignment','right');

Pulse_Dur_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string',50,'units','normalized','position',[xval 0.4 Lval 0.065],...
    'tag','Pulse_Dur');


ISI_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','ISI (msec) = ', 'units','normalized',...
    'position',[xtext 0.3 Ltext 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'horizontalalignment','right');

ISI_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string',5000,'units','normalized','position',[xval 0.3 Lval 0.065],...
    'tag','ISI');


jitter_percentage_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','Jitter (% Silence) = ', 'units','normalized',...
    'position',[xtext 0.2 Ltext 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'horizontalalignment','right');

jitter_percentage_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string',20,'units','normalized','position',[xval 0.2 Lval 0.065],...
    'tag','jitter_percentage');


Ntrials_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','Number of Trials = ', 'units','normalized',...
    'position',[xtext 0.1 Ltext 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'horizontalalignment','right');

Ntrials_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string',25,'units','normalized','position',[xval 0.1 Lval 0.065],...
    'tag','Ntrials');


Randomized_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','Randomized  ', 'units','normalized',...
    'position',[xtext 0.0 Ltext-0.17 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'ForegroundColor','b',...
    'horizontalalignment','right');

Randomized_h = uicontrol('parent',puretonesGUI,'style','checkbox',...
    'value',1,'units','normalized','position',[xval-0.175 0.002 0.1 0.08],...
    'backgroundcolor',get(puretonesGUI,'color'),'tag','Randomized');


fileno_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','File Number = ', 'units','normalized',...
    'position',[xtext+0.26 0.0 Ltext-0.16 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'horizontalalignment','right');

fileno_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string','01','units','normalized','position',[xtext+0.495 0.015 0.12 0.055],...
    'tag','fileno');


xnotes = xval + Lval + 0.05;
Notes_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','Notes: ', 'units','normalized',...
    'position',[xnotes 0.725 .1 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'ForegroundColor','b',...
    'horizontalalignment','left');


Notes_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string','','units','normalized',...
	'position',[xnotes 0.1 .4 0.65],...
	'Max',2, 'Min',0, 'horizontalalignment','left', ...
    'tag','Notes');

xSpeaker = xval + Lval + 0.05;
Speaker_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','Speaker Name:  ', 'units','normalized',...
    'position',[xSpeaker 0.9 .2 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'ForegroundColor','b',...
    'horizontalalignment','left');


Speaker_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string','ES1-A','units','normalized',...
	'position',[xSpeaker+0.15 0.91 Lval 0.055],...
	'Max',1, 'Min',0, 'horizontalalignment','center', ...
    'tag','Speaker');

xFsCam = xval + Lval + 0.05;
FsCam_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','Camera Fs (f/s): ', 'units','normalized',...
    'position',[xFsCam 0.84 .2 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'ForegroundColor','b',...
    'horizontalalignment','left');


FsCam_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string',30,'units','normalized',...
	'position',[xFsCam+0.15 0.85 Lval 0.055],...
	'Max',1, 'Min',0, 'horizontalalignment','center', ...
    'tag','FsCam');

xnFrames = xval + Lval + 0.05;
nFrames_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','nFrames: ', 'units','normalized',...
    'position',[xnFrames 0.78 .2 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'ForegroundColor','b',...
    'horizontalalignment','left');


nFrames_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string',7000,'units','normalized',...
	'position',[xnFrames+0.15 0.79 Lval 0.055],...
	'Max',1, 'Min',0, 'horizontalalignment','center', ...
    'tag','nFrames');

run_button = uicontrol('parent',puretonesGUI,'style','pushbutton',...
    'string','Run', 'units','normalized',...
    'position',[0.65 0.005 0.3 0.08],'fontsize',12,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'ForegroundColor','b',...
    'horizontalalignment','right','tag','run_button',...
    'Callback',@run_button_Callback_3tones);





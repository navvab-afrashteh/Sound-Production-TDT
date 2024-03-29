clc
clear
close all

puretonesGUI = figure('units','normalized','position',[0.2 0.5 0.45 0.45],...
    'menubar','none','numbertitle','off','name','Pure Tones Specifications');

xtext = 0.00;
Ltext = 0.39;
xval = xtext + Ltext;
Lval = 0.15;

Nvocal_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','Num Vocalization Calls = ', 'units','normalized',...
    'position',[xtext 0.9 Ltext 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'horizontalalignment','right');

Nvocal_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string',10,'units','normalized','position',[xval 0.9 Lval 0.065],...
    'tag','Nvocal');

Amp_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','Amplitude Array (dB) = ', 'units','normalized',...
    'position',[xtext 0.8 Ltext 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'horizontalalignment','right');

Amp_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string','80','units','normalized','position',[xval 0.8 1.0*Lval 0.065],...
    'fontsize',8,'tag','Amp');

ISI_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','ISI (msec) = ', 'units','normalized',...
    'position',[xtext 0.7 Ltext 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'horizontalalignment','right');

ISI_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string',3000,'units','normalized','position',[xval 0.7 Lval 0.065],...
    'tag','ISI');


jitter_percentage_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','Jitter (% Silence) = ', 'units','normalized',...
    'position',[xtext 0.6 Ltext 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'horizontalalignment','right');

jitter_percentage_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string',20,'units','normalized','position',[xval 0.6 Lval 0.065],...
    'tag','jitter_percentage');

nSponB_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','Number of Spon in the Begining = ', 'units','normalized',...
    'position',[xtext 0.5 Ltext 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'horizontalalignment','right');

nSponB_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string',2,'units','normalized','position',[xval 0.5 Lval 0.065],...
    'tag','nSponB');

nSponE_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','Number of Spon in the End = ', 'units','normalized',...
    'position',[xtext 0.4 Ltext 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'horizontalalignment','right');

nSponE_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string',0,'units','normalized','position',[xval 0.4 Lval 0.065],...
    'tag','nSponE');

Ntrials_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','Number of Trials = ', 'units','normalized',...
    'position',[xtext 0.1 Ltext 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'horizontalalignment','right');

Ntrials_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string',24,'units','normalized','position',[xval 0.1 Lval 0.065],...
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
    'position',[xnotes 0.625 .1 0.055],'fontsize',9,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'ForegroundColor','b',...
    'horizontalalignment','left');

Notes_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string','','units','normalized',...
	'position',[xnotes 0.1 .4 0.55],...
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

NCam_txt = uicontrol('parent',puretonesGUI,'style','text',...
    'string','NCam = ', 'units','normalized',...
    'position',[xnFrames 0.71 0.2 0.055],'fontsize',9,...
    'fontname','times new roman','ForegroundColor','b',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'horizontalalignment','left');

NCam_h = uicontrol('parent',puretonesGUI,'style','edit',...
    'string',5,'units','normalized',...
	'position',[xnFrames+0.15 0.72 Lval 0.055],...
	'Max',1, 'Min',0, 'horizontalalignment','center', ...
    'tag','NCam');

run_button = uicontrol('parent',puretonesGUI,'style','pushbutton',...
    'string','Run', 'units','normalized',...
    'position',[0.65 0.005 0.3 0.08],'fontsize',12,...
    'fontname','times new roman',...
    'backgroundcolor',get(puretonesGUI,'color'),...
    'ForegroundColor','b',...
    'horizontalalignment','right','tag','run_button',...
    'Callback',@run_button_Callback_Vocal_Tinnitus);






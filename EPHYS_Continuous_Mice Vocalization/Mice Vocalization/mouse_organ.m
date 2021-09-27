function [bout syllable_seq]= mouse_organ(varargin)
% function [bout syllable_seq]= mouse_organ(age_group,model_order,bout_length)
%
% Generates sequences of mouse syllabled using a first, second or third order
% Markov model to give statistically appropriate transitions between syllables
% for each age group. Sequence is output in 'bout'. In addition a waveform
% 'syllable_seq' is output constructed from typical example syllables and typical 
% intersyllable intervals for each age. These example syllables are included with
% the toolkit. A .wav file 'tempmouse.wav' of the 'syllable_seq' waveform  is 
% created in the current directory
%
% No arguments are required, the default values are:
% age_group = 0 (adults and pups)
% model_order=2; next syllable dependent on previous two syllables only
% max_bout: no bout length defined
%
% if max bout length is given, it determines the number of syllables that will
% be generated in a bout and overides the default behaviour of the model to
% end a bout only when a transition is made to an end state. This is 
% because some syllables are more likely to come at the end of a bout than 
% others.
%
% syllable types correspond to those described in Grimsley et al. 2010
% 1)	Flat syllables
% 2)	1 Frequency step syllables
% 3)	2 Frequency step syllables
% 4)	Chevron syllables
% 5)	Down-FM syllables
% 6)	Up-FM syllables
% 7)	Short syllables
% 8)	Complex syllables
% 9)	Low Frequency Harmonic syllables
% 10)	Reverse Chevron syllables
% 11)	Noisy syllables
% additional tokens are used to take account for the probabilities of each
% syllable type starting or ending a bout
% 12)   start token
% 13)   end token
%
% If age_group is given, the syllable bout produced will be statistically 
% typical for that age group
% ages groups:
% age_group = 5   : p4 pups
% age_group = 7   : p6 pups
% age_group = 9   : p8 pups
% age_group = 11  : p10 pups
% age_group = 13  : p12 pups
% age_group = 100 : adults
% age_group = 99  : all pups
% age_group = 0   : pups and adults
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%load  structure 'trans'which contains transition probability matrices for
%all age groups and model orders (e.g. trans.p5.p_1 is a 13x13 matrix of
%first order transition probabilities for p5 pups; trans.p13.p_2 a 13x13x13
%matrix of second order transition probabilities for p12 pups and so on),
%as well as matrices of number of transitions of each order (e.g.
%trans.p7.first is a 13x13 matrix of first order transitions for p7 pups)
load('transitions.mat')%transitions.mat is in the current directory

switch nargin%use default values if not given in argument
    case 0  %no arguments given use default age_group and model_order with 
            %no restriction of number of syllables in bout
        model_order=2;
        age_group = 0;
    case 1  %only age_group given, use default model_order
        model_order=2;
        age_group = varargin{1};
    case 2  %age_group and model_order given
        model_order=varargin{2};
        age_group = varargin{1};
    case 3  %age_group, model_order and bout_length given. Number of syllables 
            %in bout will be restricted to bout_length
        max_bout=varargin{3};
        model_order=varargin{2};
        age_group = varargin{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%Markov models%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch model_order
    case 1 %first order
        %initiate syllable sequence with start token:
        syllable_seq =[12];
        state=12;
        %%%%%%%%% generate rest of syllable sequence %%%%%%%%%%%%%%%%%%%%%%%%%%
        %calculate cumulative probability of transitions so random number
        % between zero and one can be used to select a transition
        cum_prob_tmp=cumsum(trans.(['p' num2str(age_group)]).p_1,2);
        if exist('max_bout','var')
            while size(syllable_seq,2)<max_bout+1%stop when desired number of syllables is reached
                %scale random number by 1-prob of end token, since the end
                %token is not used here
                if ~isnan(cum_prob_tmp(syllable_seq(end),end-1));
                a=rand*cum_prob_tmp(syllable_seq(end),end-1);
                end
                %next syllable depends only on previous syllable
                state=find(cum_prob_tmp(syllable_seq(end),:)>=a,1,'first');
                syllable_seq=[syllable_seq state];
            end
        else
            while state~=13% stop when end token is reached
                a=rand;
                %next syllable depends only on previous syllable
                state=find(cum_prob_tmp(state,:)>=a,1,'first');
                syllable_seq=[syllable_seq state];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2 % second order
        %initiate syllable sequence with start token
        syllable_seq=[12];
        %calculate cumulative probability of transitions so random number
        % between zero and one can be used to select a transition
        cum_prob_tmp=cumsum(trans.(['p' num2str(age_group)]).p_1,2);
        %need to get first syllable using first order model
        a=rand;
        state=find(cum_prob_tmp(12,:)>=a,1,'first');
        syllable_seq=[syllable_seq state];
        
        %%%%%%%%% generate rest of syllable sequence %%%%%%%%%%%%%%%%%%%%%%%%%%
        cum_prob_tmp=cumsum(trans.(['p' num2str(age_group)]).p_2,3);
        
        if exist('max_bout','var')            
            while size(syllable_seq,2)<max_bout+1%stop when desired number of syllables is reached
                %scale random number by 1-prob of end token, since the end
                %token is not used here
                if ~isnan(cum_prob_tmp(syllable_seq(end-1),syllable_seq(end),end-1))
                a=rand*cum_prob_tmp(syllable_seq(end-1),syllable_seq(end),end-1);
                end
                %next syllable depends on previous two syllables 
                state=find(cum_prob_tmp(syllable_seq(end-1),syllable_seq(end),1:end-2)>=a,1,'first');%don't allow transitions to 12 or 13
                syllable_seq=[syllable_seq state];
            end
        else
            while state~=13% stop when end token is reached
                a=rand;
                %next syllable depends on previous two calls
                state=find(cum_prob_tmp(syllable_seq(end-1),syllable_seq(end),:)>=a,1,'first');
                syllable_seq=[syllable_seq state];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3 % third order
        %initiate syllable sequence with start token
        syllable_seq=[12];
        %first syllable
        %calculate cumulative probability of transitions so random number
        %between zero and one can be used to select a transition
        cum_prob_tmp=cumsum(trans.(['p' num2str(age_group)]).p_1,2);
        a=rand;
        state=find(cum_prob_tmp(12,:)>=a,1,'first');
        syllable_seq=[syllable_seq state];
        
        %second syllable
        cum_prob_tmp=cumsum(trans.(['p' num2str(age_group)]).p_2,3);
        a=rand;
        state=find(cum_prob_tmp(12,syllable_seq(end),:)>=a,1,'first');
        syllable_seq=[syllable_seq state];
        
%       %%%%%%%%%%%%%%%%%%%%% generate rest of syllable sequence %%%%%%%%%%%%%%
        cum_prob_tmp=cumsum(trans.(['p' num2str(age_group)]).p_3,4);
       
        if exist('max_bout','var')
            while size(syllable_seq,2)<max_bout+1%stop when desired number of syllables is reached
                %scale random number by 1-prob of end token, since the end
                %token is not used here
                if ~isnan(cum_prob_tmp(syllable_seq(end-2),syllable_seq(end-1),syllable_seq(end),end-1))
                a=rand*cum_prob_tmp(syllable_seq(end-2),syllable_seq(end-1),syllable_seq(end),end-1);
                end
                %next syllable depends on previous three syllables
                state=find(cum_prob_tmp(syllable_seq(end-2),syllable_seq(end-1),syllable_seq(end),1:end-2)>=a,1,'first');
                syllable_seq=[syllable_seq state];
            end
        else
            while state~=13% stop when end token is reached
                a=rand;
                %next syllable depends on previous three syllables
                state=find(cum_prob_tmp(syllable_seq(end-2),syllable_seq(end-1),syllable_seq(end),:)>=a,1,'first');
                syllable_seq=[syllable_seq state];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Generate syllables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%concatenate typical syllables for each age group in the sequence given by the
%markov model to generate bout

%each age group has a typical inter-syllable interval in ms stored in 'gap'
%.wav files have 10 ms silence at start and end which needs to be taken
%into account
if age_group==0 ||age_group==99
    age_group=100;%use adult intervals as default
end
gap=zeros(1,100);
gap(5)=190.2-20;
gap(7)=149.1-20;
gap(9)=93.4-20;
gap(11)=96.5-20;
gap(13)=84.8-20;
gap(100)=69.4-20;

%convert sequence of syllable codes in 'syllable_seq' to syllable waveforms in 'bout'
bout=[];
%take out start and end tokens now
if exist('max_bout','var')
    syllable_seq=syllable_seq(2:end);
else
    syllable_seq=syllable_seq(2:end-1);
end

for n=1:length(syllable_seq)
        %import syllable waveform if not already done
        index=find(syllable_seq==syllable_seq(n),1,'first');
        if exist('x.index','var')
            x.n=x.index;
        else
            fname=['P' num2str(age_group) '_' num2str(syllable_seq(n)) '.wav'];
            [x.n fs]=audioread(fname);
            %insert inter-syllable intervals
            bout=[bout; zeros(round(gap(age_group)*fs/1000),1); x.n];
        end
end
%generate .wav file removing initial 10 ms silence
audiowrite('tempmouse.wav',bout(round(0.01*fs):end),fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2010 Northeastern Ohio Universities Colleges of Medicine & Pharmacy 
% and MRC Institute of Hearing Research 

% This MATLAB software was developed by Jasmine Grimsley, Jessica Monaghan
% and Jeffrey Wenstrup

% For questions regarding aspects pertaining to the vocal repertoire of mice 
% please contact Jasmine Grimsley (jgrimsley@neoucom.edu)

% For questions regarding the software please contact Jessica Monaghan
% (jessica@ihr.mrc.ac.uk)

% It is made available in the hope that it may prove useful.  Any for-
% profit use or redistribution is prohibited. No warranty is expressed
% or implied. All rights are reserved.

% Contact address:
% Dr Jasmine Grimsley,
% NEOUCOM 
% 4209 St. Rt. 44
% PO Box 95 
% Rootstown
% Ohio 44272

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



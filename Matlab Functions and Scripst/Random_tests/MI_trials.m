%% Phase-amplitude cross-frequency coupling measure
% Each event. baseline, CS-Trials, ITI

% - Performs analysis with raw and surrogate values

% The code relies on the following functions:
% --> ModIndex.m    - by Adriano Tort, Instituto do Cerebro - Universidade Federal do Rio Grande do Norte
% --> shuffle_esc.m - by Rafal Bogacz, Angela Onslow, May 2010

% See Tort et al, 2010 -> 10.1152/jn.00106.2010 
%     Tort et al, 2008 -> 10.1073/pnas.0810524105

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  01/2024
% Last update: 01/2024

%% Hilbert Transform

% MI.data

% 1st column - > phase
% 2nd column - > amplitude
% 3nd column - > time vector

% in each cell -> frequencies cutoff
% - rows          - > channels
% - columns       - > time
% - 3nd dimension - > behavioral events

% Choose frequency bands index according to the Pre_processing.m
% MI.params.freq_idx:
% 2 - longcutoff / 3 - theta / 4 - low heta (3-6)
% 5  - hightheta (6 - 9) / 6 - lowgamma (30 50)  / 7 - highgamma (62 100)
% 8  - theta (5 7) / 9 = gamma (80 90)

MI.params.freq_idx = [8 9]; %frequency 

% Loop over channels and make hilbert transform
% Initializing with NaN
MI.data = cell(3,size(MI.params.freq_idx,2));

MI.data(1,:) = {nan(size(data.lfp{6,1}))}; % baseline
MI.data(2,:) = {nan(size(data.lfp{7,1}))}; % CS-Trials
MI.data(3,:) = {nan(size(data.lfp{8,1}))}; % ITI

% baseline
for ii = 1:length(MI.params.freq_idx)
    for jj = 1:size(data.lfp{6, MI.params.freq_idx(ii)},3)
        for ll = 1:size(data.lfp{6, MI.params.freq_idx(ii)},1)
        
            temp = data.lfp{6, MI.params.freq_idx(ii)}(ll,:,jj);
            temp(isnan(temp(:,:)))=[]; % it was necessary to remove the nan for the function to work

            if ii == 1
                MI.data{1, ii}(ll,1:length(temp),jj) = angle(hilbert(temp));
            else
                MI.data{1, ii}(ll,1:length(temp),jj) = abs(hilbert(temp));
            end
    
        end 
    end 
end


% CS-trials
for ii = 1:length(MI.params.freq_idx)
    for jj = 1:size(data.lfp{7, MI.params.freq_idx(ii)},3)
        for ll = 1:size(data.lfp{7, MI.params.freq_idx(ii)},1)
        
            temp = data.lfp{7, MI.params.freq_idx(ii)}(ll,:,jj);
            temp(isnan(temp(:,:)))=[]; % it was necessary to remove the nan for the function to work

            if ii == 1
                MI.data{2, ii}(ll,1:length(temp),jj) = angle(hilbert(temp));
            else
                MI.data{2, ii}(ll,1:length(temp),jj) = abs(hilbert(temp));
            end
    
        end 
    end 
end


% ITI
for ii = 1:length(MI.params.freq_idx)
    for jj = 1:size(data.lfp{8, MI.params.freq_idx(ii)},3)
        for ll = 1:size(data.lfp{8, MI.params.freq_idx(ii)},1)
        
            temp = data.lfp{8, MI.params.freq_idx(ii)}(ll,:,jj);
            temp(isnan(temp(:,:)))=[]; % it was necessary to remove the nan for the function to work

            if ii == 1
                MI.data{3, ii}(ll,1:length(temp),jj) = angle(hilbert(temp));
            else
                MI.data{3, ii}(ll,1:length(temp),jj) = abs(hilbert(temp));
            end
    
        end 
    end 
end


% Time vector
MI.data{1,3} = linspace(0,size(MI.data{1,1},2)/parameters.decimated_srate,size(MI.data{1, 1},2));
MI.data{2,3} = linspace(0,size(MI.data{2,1},2)/parameters.decimated_srate,size(MI.data{2, 1},2));
MI.data{3,3} = linspace(0,size(MI.data{3,1},2)/parameters.decimated_srate,size(MI.data{3, 1},2));


clear('ii','jj','ll','temp')

%% Modulation Index

% ----------------------------------------------
% IMPORTANTE NOTE:
%
% The phase and amplitude values for the MI calculation are being taken from the same channel!!!!!!!!
% The code needs to be modified for multiple comparisons between channels.
% ----------------------------------------------


% MI_values 
% - rows          - > 1st: baseline
%                 - > 2st: CS-trials
%                 - > 3st: ITI
% -each cell      - > Lines: channels
%                 - > Columns: MI values


% MeanAmp 
% - rows          - > 1st: baseline
%                 - > 2st: CS-trials
%                 - > 3st: ITI
% - each cell     - > Lines: channels
%                 - > columns:  amplitude distribution over phase bins


% Define time to analyse -> not prepare for this...
% MI.params.time2analise = 2;
% time2samples = MI.params.time2analise * parameters.decimated_srate; % in samples


% Define number number of phase bins 
MI.params.nbin = 40; 

% variable not centered in the phase bin (rad)
MI.params.position = zeros(1,MI.params.nbin); 

MI.params.winsize = 2*pi/MI.params.nbin;

for jj = 1:MI.params.nbin
    MI.params.position(jj) = -pi+(jj-1)*MI.params.winsize;
end


% Baseline
MI.MI_value{1,1} = nan(size(MI.data{1,1},1),1,size(MI.data{1,1},3));
MI.MeanAmp{1,1}  = nan(size(MI.data{1,1},1),MI.params.nbin,size(MI.data{1,1},3));

for ii = 1:size(MI.data{1,1},1)
    [MI.MI_value{1,1}(ii,1),MI.MeanAmp{1,1}(ii,:)] = ModIndex(MI.data{1,1}(ii,:), MI.data{1,2}(ii,:), MI.params.position);
end


% CS-Trial
MI.MI_value{2,1} = nan(size(MI.data{2,1},1),size(MI.data{2,1},3));
MI.MeanAmp{2,1}  = nan(size(MI.data{2,1},1),MI.params.nbin,size(MI.data{2,1},3));

for ii = 1:size(MI.data{2,1},1)
    for jj = 1:size(MI.data{2,1},3)

        [MI.MI_value{2,1}(ii,jj),MI.MeanAmp{2,1}(ii,:,jj)] = ModIndex(MI.data{2,1}(ii,:,jj), MI.data{2,2}(ii,:,jj), MI.params.position);

    end
end


% ITI
MI.MI_value{3,1} = nan(size(MI.data{3,1},1),size(MI.data{3,1},3));
MI.MeanAmp{3,1}  = nan(size(MI.data{3,1},1),MI.params.nbin,size(MI.data{3,1},3));

for ii = 1:size(MI.data{3,1},1)
    for jj = 1:size(MI.data{3,1},3)

        [MI.MI_value{3,1}(ii,jj),MI.MeanAmp{3,1}(ii,:,jj)] = ModIndex(MI.data{3,1}(ii,:,jj), MI.data{3,2}(ii,:,jj), MI.params.position);

    end
end


clear('ii','jj','time2samples')

%% Plot MI values

%Choose channel: 
% ch_phase = 6
% ch_amp   = 6

ch = 6;

figure
set(gcf,'color','w');

% sgtitle({' ';'\fontsize{14}Phase-amplitude cross-frequency coupling'; 
%     ['Phase channel: ' num2str(ch_phase),' / Amplitude channel: ' num2str(ch_amp)];' ' });

sgtitle('\fontsize{14}Phase-amplitude cross-frequency coupling');

subplot 221
b1 = bar(MI.MI_value{1,1}(ch,:));
b1(1).FaceColor = 'w';
ylim([0 1.1*10^-3])

for ii = 1:size(MI.MI_value{1,1},2)
    labels(ii)={['Baseline: ' num2str(ii)]};
end

set(gca,'xticklabel',labels(1:1:end));
ylabel('\fontsize{14} Modulation Index');


subplot 222
b1 = bar(MI.MI_value{2,1}(ch,:));
b1(1).FaceColor = 'w';
ylim([0 1.1*10^-3])

for ii = 1:size(MI.MI_value{2,1},2)
    labels(ii)={['CS: ' num2str(ii)]};
end

set(gca,'xticklabel',labels(1:1:end));
ylabel('\fontsize{14} Modulation Index');


subplot 224
b3 = bar(MI.MI_value{3,1}(ch,:));
b3(1).FaceColor = 'w';
ylim([0 1.1*10^-3])

for ii = 1:size(MI.MI_value{3,1},2)
    labels(ii)={['ITI: ' num2str(ii)]};
end

set(gca,'xticklabel',labels(1:1:end));
ylabel('\fontsize{14} Modulation Index');



clear('b1','b2','b3','labels','ii','ch')

%% Plot phase x amplitude distributions

ch = 6;

figure
set(gcf,'color','w');

xvalue1 = (rad2deg(MI.params.position));% + size(MI.params.position,2)/2);

    
subplot(2,size(MI.MeanAmp{2,1},3)+1,1)
b3 = bar([xvalue1] ,[MI.MeanAmp{1,1}(ch,:)]);
b3.FaceColor = 'w';
ylabel('\fontsize{11} amplitude');
xlabel('\fontsize{11} phase');
box off
title('baseline')


for ii = 1: size(MI.MeanAmp{2,1},3)
    subplot(2,size(MI.MeanAmp{2,1},3)+1,ii+1)
    b2 = bar(xvalue1,MI.MeanAmp{2,1}(ch,:,ii));
    b2.FaceColor = 'w';
    ylabel('\fontsize{11} amplitude');
    xlabel('\fontsize{11} phase');
    box off

    title ({['CS: ' num2str(ii)]})

end
  

for ii = 1: size(MI.MeanAmp{3,1},3)
    subplot(2,size(MI.MeanAmp{3,1},3)+1,ii+7)
    b3 = bar(xvalue1,MI.MeanAmp{3,1}(ch,:,ii));
    b3.FaceColor = 'w';
    ylabel('\fontsize{11} amplitude');
    xlabel('\fontsize{11} phase');
    box off

    title ({['ITI: ' num2str(ii)]})

end
    
clear('xvalue1','ii','b2','b3','b1','ch')

%% Plot Theta and Gamma Oscillation

ch = 6;

% Time vector
timev_baseline  = linspace(0,length(data.lfp{6,1}(ch,:))./parameters.decimated_srate,length(data.lfp{6,1}(ch,:)));
timev_CS_Trials = linspace(0,length(data.lfp{7,1}(ch,:,1))./parameters.decimated_srate,length(data.lfp{7,1}(ch,:,1)));
timev_ITI       = linspace(0,length(data.lfp{8,1}(ch,:,1))./parameters.decimated_srate,length(data.lfp{8,1}(ch,:,1)));


figure
set(gcf,'color','w');
    
subplot(2,size(MI.MeanAmp{2,1},3)+1,1)

plot(timev_baseline,data.lfp{6,MI.params.freq_idx(1,2)}(ch,:),'Color',[.6, 0, 0 .7])
hold on
plot(timev_baseline,data.lfp{6,MI.params.freq_idx(1,1)}(ch,:),'linew',1.5, 'Color',[.2, .2, .2])

ylabel('\fontsize{11} uV');
xlabel('\fontsize{11} Time (s)');
xlim([0 timev_baseline(end)])
box off
title('baseline')


for ii = 1: size(MI.MeanAmp{2,1},3)
    subplot(2,size(MI.MeanAmp{2,1},3)+1,ii+1)
    plot(timev_CS_Trials,data.lfp{7,MI.params.freq_idx(1,2)}(ch,:,ii),'Color',[.6, 0, 0 .7])
    hold on
    plot(timev_CS_Trials,data.lfp{7,MI.params.freq_idx(1,1)}(ch,:,ii),'linew',1.5, 'Color',[.2, .2, .2])

    ylabel('\fontsize{11} uV');
    xlabel('\fontsize{11} Time (s)');
    xlim([0 timev_CS_Trials(end)])
    box off

    title ({['CS: ' num2str(ii)]})

end
  

for ii = 1: size(MI.MeanAmp{3,1},3)
    subplot(2,size(MI.MeanAmp{2,1},3)+1,ii+7)
    plot(timev_ITI,data.lfp{8,MI.params.freq_idx(1,2)}(ch,:,ii),'Color',[.6, 0, 0 .7])
    hold on
    plot(timev_ITI,data.lfp{8,MI.params.freq_idx(1,1)}(ch,:,ii),'linew',1.5, 'Color',[.2, .2, .2])

    ylabel('\fontsize{11} uV');
    xlabel('\fontsize{11} Time (s)');
    xlim([0 timev_ITI(end)])
    box off

    title ({['ITI: ' num2str(ii)]})

end
    
clear('ch','ii','timev_baseline','timev_CS_Trials','timev_ITI')


%% Surrogate phase vectors to compare

% MI_shuffle_values 
% - rows          - > events
% - 1st column    - > before onset
% - 2nd column    - > post onset
% - 3nd dimension - > rearrangements

% MeanAmp_shuffle 
% - rows          - > events
% - coluns        - > amplitude distribution over phase bins
% - 3nd dimension - > 1st: before onset
%                     2nd: post onset
% - 4nd dimension - > rearrangements

numshf  = 200; % number of shuffled segments
nsurrog = 200; % number of rearrangements


% Loop over behavior events 
MI.MI_shuffle_values = [];
MI.MeanAmp_shuffle   = [];

% baseline
baseline_shuffle = [];
for jj = 1:nsurrog
    
    for ii = 1:size(MI.data{1,1},3)
    
        baseline_shuffle(:,:,ii) = shuffle_esc(MI.data{1,1}(:,:,ii),parameters.decimated_srate,numshf);
    
        for cc = 1:size(MI.data{1,1},1)

            [MI.MI_shuffle_values{1,1}(cc,jj,ii),MI.MeanAmp_shuffle{1,1}(cc,:,jj,ii)] = ModIndex(baseline_shuffle(cc,:,ii), MI.data{1,2}(cc,:,ii), MI.params.position);

        end    
    end    
end

% CS-Trials
CS_shuffle = [];
for jj = 1:nsurrog
    
    for ii = 1:size(MI.data{2,1},3)
    
        CS_shuffle(:,:,ii) = shuffle_esc(MI.data{2,1}(:,:,ii),parameters.decimated_srate,numshf);
    
        for cc = 1:size(MI.data{2,1},1)

            [MI.MI_shuffle_values{2,1}(cc,jj,ii),MI.MeanAmp_shuffle{2,1}(cc,:,jj,ii)] = ModIndex(CS_shuffle(cc,:,ii), MI.data{2,2}(cc,:,ii), MI.params.position);

        end    
    end    
end

% ITI
ITI_shuffle = [];
for jj = 1:nsurrog
    
    for ii = 1:size(MI.data{3,1},3)
    
        ITI_shuffle(:,:,ii) = shuffle_esc(MI.data{3,1}(:,:,ii),parameters.decimated_srate,numshf);
    
        for cc = 1:size(MI.data{3,1},1)

            [MI.MI_shuffle_values{3,1}(cc,jj,ii),MI.MeanAmp_shuffle{3,1}(cc,:,jj,ii)] = ModIndex(ITI_shuffle(cc,:,ii), MI.data{3,2}(cc,:,ii), MI.params.position);

        end    
    end    
end

% z surrogated values
MI.stats.z_MI_shuffle_values{1,1} = zscore(MI.MI_shuffle_values{1,1});
MI.stats.z_MI_shuffle_values{2,1} = zscore(MI.MI_shuffle_values{2,1});
MI.stats.z_MI_shuffle_values{3,1} = zscore(MI.MI_shuffle_values{3,1});


% z real(observed) values
MI.stats.z_MI_value{1,1} = (MI.MI_value{1,1} - mean(MI.MI_shuffle_values{1,1},2))./squeeze(std(MI.MI_shuffle_values{1,1},[],2));
MI.stats.z_MI_value{2,1} = (MI.MI_value{2,1} - squeeze(mean(MI.MI_shuffle_values{2,1},2)))./squeeze(std(MI.MI_shuffle_values{2,1},[],2));
MI.stats.z_MI_value{3,1} = (MI.MI_value{3,1} - squeeze(mean(MI.MI_shuffle_values{3,1},2)))./squeeze(std(MI.MI_shuffle_values{3,1},[],2));

% p values
MI.stats.p_MI_value{1,1} = normcdf(-MI.stats.z_MI_value{1,1});
MI.stats.p_MI_value{2,1} = normcdf(-MI.stats.z_MI_value{2,1});
MI.stats.p_MI_value{3,1} = normcdf(-MI.stats.z_MI_value{3,1});

%MI.stats.p_MI_value = 0.5 * erfc(MI.stats.z_MI_value ./ sqrt(2));
%MI.stats.p_MI_value = erfc(- MI.stats.z_MI_value ./ sqrt(2))./2;


clear('jj','ii','cc','nsurrog','numshf','time2samples','baseline_shuffle','CS_shuffle','ITI_shuffle')

%% Plot statistical distribuitions

ch = 6;

figure
set(gcf,'color','w');

MI.params.nbins = 30;



% Baseline
subplot(2,size(MI.MeanAmp{2,1},3)+1,1)
h = histogram(MI.stats.z_MI_shuffle_values{1,1}(ch,:),MI.params.nbins);
h.FaceColor = 'w';
hold on
plot([MI.stats.z_MI_value{1,1}(ch,:) MI.stats.z_MI_value{1,1}(ch,:)],[0 30],'Color', '[0.6350, 0.0780, 0.1840]','linewidth',4);
xlim([-2.5 2.5])
%ylim ([0 30])

% just to make arrow
%[figx,figy] = dsxy2figxy([MI.stats.z_MI_value(ii,1) MI.stats.z_MI_value(ii,1)],[10 0]); % Transform point or MI.params.position from data space
%                                                                                          coordinates into normalized figure coordinates .
%                                                                                          >>> https://uk.mathworks.com/matlabcentral/fileexchange/30347-sigplot?focused=5178148&tab=function
%annotation('textarrow',figx,figy,'Color', '[0.6350, 0.0780, 0.1840]','linewidth',4);

title({'Habituation: ';['MI = ' num2str(MI.MI_value{1,1}(ch,1)*10^3),'  ','\rm\fontsize{12}p_v_a_l_u_e = ' num2str(MI.stats.p_MI_value{1,1}(ch,1))]})

ylabel('\fontsize{14}Frequency')
xlabel('\fontsize{12}z values')
%title(['MI = ' num2str(MI.MI_value(count,1)),'  ','\rm\fontsize{12}p_v_a_l_u_e = ' num2str(MI.stats.p_MI_value(count,1)) ])
%legend({'\fontsize{12}Random';'\fontsize{12}Observed'},'box','off');

box off

%CS-Trials
for ii = 1: size(MI.MeanAmp{2,1},3)
    subplot(2,size(MI.MeanAmp{2,1},3)+1,ii+1)
    
    h = histogram(MI.stats.z_MI_shuffle_values{2,1}(ch,:,ii),MI.params.nbins);
    h.FaceColor = 'w';
    hold on
    plot([MI.stats.z_MI_value{2,1}(ch,ii) MI.stats.z_MI_value{2,1}(ch,ii)],[0 30],'Color', '[0.6350, 0.0780, 0.1840]','linewidth',4);
    xlim([-2.5 2.5])
    %ylim ([0 30])
    
    % just to make arrow
    %[figx,figy] = dsxy2figxy([MI.stats.z_MI_value(ii,1) MI.stats.z_MI_value(ii,1)],[10 0]); % Transform point or MI.params.position from data space 
                                                                                            % coordinates into normalized figure coordinates . 
                                                                                            % >>> https://uk.mathworks.com/matlabcentral/fileexchange/30347-sigplot?focused=5178148&tab=function
    %annotation('textarrow',figx,figy,'Color', '[0.6350, 0.0780, 0.1840]','linewidth',4);
    
    title({['CS-Trials: ' num2str(ii)];['MI = ' num2str(MI.MI_value{2,1}(ch,ii)*10^3),'  ','\rm\fontsize{12}p_v_a_l_u_e = ' num2str(MI.stats.p_MI_value{2,1}(ch,ii))]})    
    
    ylabel('\fontsize{14}Frequency')
    xlabel('\fontsize{12}z values')
    %title(['MI = ' num2str(MI.MI_value(count,1)),'  ','\rm\fontsize{12}p_v_a_l_u_e = ' num2str(MI.stats.p_MI_value(count,1)) ])
    %legend({'\fontsize{12}Random';'\fontsize{12}Observed'},'box','off');

    box off

end


%ITI
for ii = 1: size(MI.MeanAmp{3,1},3)
    subplot(2,size(MI.MeanAmp{3,1},3)+1,ii+7)
    
    h = histogram(MI.stats.z_MI_shuffle_values{3,1}(ch,:,ii),MI.params.nbins);
    h.FaceColor = 'w';
    hold on
    plot([MI.stats.z_MI_value{3,1}(ch,ii) MI.stats.z_MI_value{3,1}(ch,ii)],[0 30],'Color', '[0.6350, 0.0780, 0.1840]','linewidth',4);
    xlim([-2.5 2.5])
    %ylim ([0 30])
    
    % just to make arrow
    %[figx,figy] = dsxy2figxy([MI.stats.z_MI_value(ii,1) MI.stats.z_MI_value(ii,1)],[10 0]); % Transform point or MI.params.position from data space 
                                                                                            % coordinates into normalized figure coordinates . 
                                                                                            % >>> https://uk.mathworks.com/matlabcentral/fileexchange/30347-sigplot?focused=5178148&tab=function
    %annotation('textarrow',figx,figy,'Color', '[0.6350, 0.0780, 0.1840]','linewidth',4);
    

    title({['ITI: ' num2str(ii)];['MI = ' num2str(MI.MI_value{3,1}(ch,ii)*10^3),'  ','\rm\fontsize{12}p_v_a_l_u_e = ' num2str(MI.stats.p_MI_value{3,1}(ch,ii))]})    
    
    ylabel('\fontsize{14}Frequency')
    xlabel('\fontsize{12}z values')
    %title(['MI = ' num2str(MI.MI_value(count,1)),'  ','\rm\fontsize{12}p_v_a_l_u_e = ' num2str(MI.stats.p_MI_value(count,1)) ])
    %legend({'\fontsize{12}Random';'\fontsize{12}Observed'},'box','off');

    box off

end

clear ('ii','jj','count','h','figx','figy','ch');

%% last update 14/01/2024 - 18:51
%  listening: PJ harvey - White Chalk

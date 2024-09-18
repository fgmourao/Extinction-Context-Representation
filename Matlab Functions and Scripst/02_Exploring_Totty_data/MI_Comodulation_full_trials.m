%% Phase-amplitude Cross-frequency coupling measure

% Based on the Hindiael Belchior & Adriano Tort script
% Instituto do Cerebro - Universidade Federal do Rio Grande do Norte.Brazil


% The code relies on the following functions:
% --> ModIndex.m - by Adriano Tort, Instituto do Cerebro - Universidade Federal do Rio Grande do Norte

% See Tort et al, 2010 -> 10.1152/jn.00106.2010 
%     Tort et al, 2008 -> 10.1073/pnas.0810524105


% Adapted by:
% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  01/2024
% Last update: 04/2024

%%
% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m

%% Run each session sequentially

MI.comod = [];

% Pair channels combination
% Row 1: PL  --> IL
% Row 2: IL  --> PL
% Row 3: PL  --> HPC
% Row 4: HPC --> PL
% Row 5: IL  --> HPC
% Row 6: HPC --> IL

% All Possibles combinations
combinations_ = zeros(size(data.lfp{6,1},1)*3,2);
combinations_(1:2:6,:) = flip(nchoosek(1:size(data.lfp{6,1},1),2),2);
combinations_(2:2:6,:) = nchoosek(1:size(data.lfp{6,1},1),2);
combinations_(7:9,1) = 1:size(data.lfp{6,1},1);
combinations_(7:9,2) = 1:size(data.lfp{6,1},1);

% Choose channels
MI.comod.parameters.ch_phase = combinations_(:,1);
MI.comod.parameters.ch_amp   = combinations_(:,2);

for ii = 1:size(combinations_,1)
    MI.comod.lfp_phase_baseline{ii,1} = data.lfp{6,1}(MI.comod.parameters.ch_phase(ii,1),B_clean{ms}(1):B_clean{ms}(2));
    MI.comod.lfp_amp_baseline{ii,1}   = data.lfp{6,1}(MI.comod.parameters.ch_amp(ii,1),B_clean{ms}(1):B_clean{ms}(2));

    MI.comod.lfp_phase_CS{ii,1}       = data.lfp{7,1}(MI.comod.parameters.ch_phase(ii,1),:,CSIT{ms});
    MI.comod.lfp_amp_CS{ii,1}         = data.lfp{7,1}(MI.comod.parameters.ch_amp(ii,1),:,CSIT{ms});

    MI.comod.lfp_phase_ITI{ii,1}      = data.lfp{8,1}(MI.comod.parameters.ch_phase(ii,1),:,CSIT_1{ms});
    MI.comod.lfp_amp_ITI{ii,1}        = data.lfp{8,1}(MI.comod.parameters.ch_amp(ii,1),:,CSIT_1{ms});
end

clear('ii')

%% Setings - Define the Amplitude and Phase Frequencies range

% MI.comod.parameters.PhaseFreqVector = 1:.5:20;
% MI.comod.parameters.AmpFreqVector   = 30:2:120;
% 
% MI.comod.parameters.PhaseFreq_BandWidth = 1;
% MI.comod.parameters.AmpFreq_BandWidth   = 5;

%Default by author
MI.comod.parameters.PhaseFreqVector = 2.5:.5:12;
MI.comod.parameters.AmpFreqVector   = 20:2:150;

MI.comod.parameters.PhaseFreq_BandWidth = 2;
MI.comod.parameters.AmpFreq_BandWidth   = 5;


%% Do filtering and Hilbert transform on CPU

'CPU filtering';

tic
disp('Filtering and Hilbert transform loop on baseline')

for cc = 1:size(combinations_,1)

    % baseline

    MI.comod.AmpFreqTransformed_baseline{cc,1}   = zeros(length(MI.comod.parameters.AmpFreqVector), length(MI.comod.lfp_amp_baseline{cc,1}),size(MI.comod.lfp_amp_baseline{cc,1},3));
    MI.comod.PhaseFreqTransformed_baseline{cc,1} = zeros(length(MI.comod.parameters.PhaseFreqVector), length(MI.comod.lfp_phase_baseline{cc,1}), size(MI.comod.lfp_phase_baseline{cc,1},3));

    for jj = 1:size(MI.comod.lfp_amp_baseline{cc,1} ,3)

        for ii=1:length(MI.comod.parameters.AmpFreqVector)
            Af1 = MI.comod.parameters.AmpFreqVector(ii);
            Af2 = Af1+MI.comod.parameters.AmpFreq_BandWidth;
            AmpFreq = eegfilt(MI.comod.lfp_amp_baseline{cc,1}(:,:,jj),parameters.decimated_srate,Af1,Af2); % just filtering
            MI.comod.AmpFreqTransformed_baseline{cc,1}(ii,:,jj) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
        end

    end


    for jj = 1:size(MI.comod.lfp_phase_baseline{cc,1} ,3)
        for ii=1:length(MI.comod.parameters.PhaseFreqVector)
            Pf1 = MI.comod.parameters.PhaseFreqVector(ii);
            Pf2 = Pf1 + MI.comod.parameters.PhaseFreq_BandWidth;
            PhaseFreq = eegfilt(MI.comod.lfp_phase_baseline{cc,1}(:,:,jj),parameters.decimated_srate,Pf1,Pf2); % this is just filtering
            MI.comod.PhaseFreqTransformed_baseline{cc,1}(ii,:,jj) = angle(hilbert(PhaseFreq)); % this is getting the phase time series
        end
    end


end

toc


tic
disp('Filtering and Hilbert transform loop on CS-Trials')

for cc = 1:size(combinations_,1)

    % CS-Trials

    MI.comod.AmpFreqTransformed_CS{cc,1}    = zeros(length(MI.comod.parameters.AmpFreqVector), length(MI.comod.lfp_amp_CS{cc,1}),size(MI.comod.lfp_amp_CS{cc,1},3));
    MI.comod.PhaseFreqTransformed_CS{cc,1}  = zeros(length(MI.comod.parameters.PhaseFreqVector), length(MI.comod.lfp_phase_CS{cc,1}), size(MI.comod.lfp_phase_CS{cc,1},3));

    for jj = 1:size(MI.comod.lfp_amp_CS{cc,1}  ,3)
        for ii=1:length(MI.comod.parameters.AmpFreqVector)
            Af1 = MI.comod.parameters.AmpFreqVector(ii);
            Af2 = Af1+MI.comod.parameters.AmpFreq_BandWidth;
            AmpFreq = eegfilt(MI.comod.lfp_amp_CS{cc,1}(:,:,jj),parameters.decimated_srate,Af1,Af2); % just filtering
            MI.comod.AmpFreqTransformed_CS{cc,1}(ii,:,jj) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
        end
    end

    for jj = 1:size(MI.comod.lfp_phase_CS{cc,1} ,3)
        for ii=1:length(MI.comod.parameters.PhaseFreqVector)
            Pf1 = MI.comod.parameters.PhaseFreqVector(ii);
            Pf2 = Pf1 + MI.comod.parameters.PhaseFreq_BandWidth;
            PhaseFreq = eegfilt(MI.comod.lfp_phase_CS{cc,1}(:,:,jj),parameters.decimated_srate,Pf1,Pf2); % this is just filtering
            MI.comod.PhaseFreqTransformed_CS{cc,1}(ii,:,jj) = angle(hilbert(PhaseFreq)); % this is getting the phase time series
        end
    end

end

toc


tic
disp('Filtering and Hilbert transform loop on ITIs')

for cc = 1:size(combinations_,1)

    % ITI

    MI.comod.AmpFreqTransformed_ITI{cc,1}   = zeros(length(MI.comod.parameters.AmpFreqVector), length(MI.comod.lfp_amp_ITI{cc,1}),size(MI.comod.lfp_amp_ITI{cc,1},3));
    MI.comod.PhaseFreqTransformed_ITI{cc,1} = zeros(length(MI.comod.parameters.PhaseFreqVector), length(MI.comod.lfp_phase_ITI{cc,1}), size(MI.comod.lfp_phase_ITI{cc,1},3));

    for jj = 1:size(MI.comod.lfp_amp_ITI{cc,1} ,3)
        for ii=1:length(MI.comod.parameters.AmpFreqVector)
            Af1 = MI.comod.parameters.AmpFreqVector(ii);
            Af2 = Af1+MI.comod.parameters.AmpFreq_BandWidth;
            AmpFreq = eegfilt(MI.comod.lfp_amp_ITI{cc,1}(:,:,jj),parameters.decimated_srate,Af1,Af2); % just filtering
            MI.comod.AmpFreqTransformed_ITI{cc,1}(ii,:,jj) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
        end
    end

    for jj = 1:size(MI.comod.lfp_phase_ITI{cc,1} ,3)
        for ii=1:length(MI.comod.parameters.PhaseFreqVector)
            Pf1 = MI.comod.parameters.PhaseFreqVector(ii);
            Pf2 = Pf1 + MI.comod.parameters.PhaseFreq_BandWidth;
            PhaseFreq = eegfilt(MI.comod.lfp_phase_ITI{cc,1}(:,:,jj),parameters.decimated_srate,Pf1,Pf2); % this is just filtering
            MI.comod.PhaseFreqTransformed_ITI{cc,1}(ii,:,jj) = angle(hilbert(PhaseFreq)); % this is getting the phase time series
        end
    end


end

toc

clear('Trial_','ans','cc','ii', 'jj', 'Af1', 'Af2', 'Pf1', 'Pf2', 'AmpFreq', 'PhaseFreq');


% For comodulation calculation (only has to be calculated once)

MI.comod.parameters.nbin     = 18;
MI.comod.parameters.position = zeros(1,MI.comod.parameters.nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
MI.comod.parameters.winsize  = 2*pi/MI.comod.parameters.nbin;

for j=1:MI.comod.parameters.nbin 
    MI.comod.parameters.position(j) = -pi+(j-1)*MI.comod.parameters.winsize; 
end


% Comodulogram
% Cell rows -> channels combinations

% each cell
%   - rows            -> phase
%   - columns         -> amplitude
%   - third dimension -> events


% baseline
tic
disp('Comodulation loop baseline')

for cc = 1:size(combinations_,1)

    MI.comod.Comodulogram_baseline{cc,1} = single(zeros(length(MI.comod.parameters.PhaseFreqVector),length(MI.comod.parameters.AmpFreqVector)));

    for nn = 1:size(MI.comod.PhaseFreqTransformed_baseline{cc,1},3)
        
        counter1=0;
        for ii = 1:length(MI.comod.parameters.PhaseFreqVector)
            counter1 = counter1 + 1;

            MI.comod.Pf1_baseline(ii) = MI.comod.parameters.PhaseFreqVector(ii);
            MI.comod.Pf2_baseline     = MI.comod.Pf1_baseline + MI.comod.parameters.PhaseFreq_BandWidth;

            counter2=0;
            for jj = 1:length(MI.comod.parameters.AmpFreqVector)
                counter2 = counter2 + 1;

                MI.comod.Af1_baseline(jj) = MI.comod.parameters.AmpFreqVector(jj);
                MI.comod.Af2_baseline     = MI.comod.Af1_baseline + MI.comod.parameters.AmpFreq_BandWidth;

                [MI_baseline,~] = ModIndex(MI.comod.PhaseFreqTransformed_baseline{cc,1}(ii,:,nn), MI.comod.AmpFreqTransformed_baseline{cc,1}(jj,:,nn), MI.comod.parameters.position);

                MI.comod.Comodulogram_baseline{cc,1}(counter1,counter2,nn) = MI_baseline;

            end
        end
    end
end

toc


% CS-Trials
tic
disp('Comodulation loop CS-Trials')

for cc = 1:size(combinations_,1)

    MI.comod.Comodulogram_CS{cc,1} = single(zeros(length(MI.comod.parameters.PhaseFreqVector),length(MI.comod.parameters.AmpFreqVector),CS_));

    for nn = 1:size(MI.comod.PhaseFreqTransformed_CS{cc,1},3)
        
        counter1=0;
        for ii = 1:length(MI.comod.parameters.PhaseFreqVector)
            counter1 = counter1 + 1;

            MI.comod.Pf1_CS(ii) = MI.comod.parameters.PhaseFreqVector(ii);
            MI.comod.Pf2_CS     = MI.comod.Pf1_CS + MI.comod.parameters.PhaseFreq_BandWidth;

            counter2=0;
            for jj = 1:length(MI.comod.parameters.AmpFreqVector)
                counter2 = counter2 + 1;

                MI.comod.Af1_CS(jj) = MI.comod.parameters.AmpFreqVector(jj);
                MI.comod.Af2_CS     = MI.comod.Af1_CS + MI.comod.parameters.AmpFreq_BandWidth;

                [MI_CS,~] = ModIndex(MI.comod.PhaseFreqTransformed_CS{cc,1}(ii,:,nn), MI.comod.AmpFreqTransformed_CS{cc,1}(jj,:,nn), MI.comod.parameters.position);

                MI.comod.Comodulogram_CS{cc,1}(counter1,counter2,nn) = MI_CS;

            end
        end
    end
end

toc


% ITI
tic
disp('Comodulation loop ITIs')

for cc = 1:size(combinations_,1)

    MI.comod.Comodulogram_ITI{cc,1} = single(zeros(length(MI.comod.parameters.PhaseFreqVector),length(MI.comod.parameters.AmpFreqVector),ITI_));

    for nn = 1:size(MI.comod.PhaseFreqTransformed_ITI{cc,1},3)
        
        counter1=0;
        for ii = 1:length(MI.comod.parameters.PhaseFreqVector)
            counter1 = counter1 + 1;

            MI.comod.Pf1_ITI(ii) = MI.comod.parameters.PhaseFreqVector(ii);
            MI.comod.Pf2_ITI     = MI.comod.Pf1_ITI + MI.comod.parameters.PhaseFreq_BandWidth;

            counter2=0;
            for jj = 1:length(MI.comod.parameters.AmpFreqVector)
                counter2 = counter2 + 1;

                MI.comod.Af1_ITI(jj) = MI.comod.parameters.AmpFreqVector(jj);
                MI.comod.Af2_ITI     = MI.comod.Af1_ITI + MI.comod.parameters.AmpFreq_BandWidth;

                [MI_ITI,~] = ModIndex(MI.comod.PhaseFreqTransformed_ITI{cc,1}(ii,:,nn), MI.comod.AmpFreqTransformed_ITI{cc,1}(jj,:,nn), MI.comod.parameters.position);

                MI.comod.Comodulogram_ITI{cc,1}(counter1,counter2,nn) = MI_ITI;

            end
        end
    end
end

toc

disp('Done')

clear ('cc','j','MI_CS','MI_ITI','ii','jj','nn','time2samples','MI_baseline','counter1','counter2')

%% Graph comodulogram
%  Events pars, avereged trials - Baseline, CS-Tones and ITI


toplot_       = cell(max(max(combinations_)),size(combinations_,1));
toplot_(1,:)  = MI.comod.Comodulogram_baseline;
toplot_(2,:)  = MI.comod.Comodulogram_CS;
toplot_(3,:)  = MI.comod.Comodulogram_ITI;

titles_ = { 'Baseline','CS-Tones','ITI'};
y_combination = {'IL phase \nPL ampl','PL phase - IL ampl',...
    'dHPC phase - PL ampl', 'PL phase - dHPC ampl',...
    'dHPC phase - IL ampl', 'IL phase - dHPC ampl',...
    'PL phase - PL ampl',...
    'IL phase - IL ampl',...
    'dHPC phase - dHPC ampl'};


figure
set(gcf,'color','w');
sc = [1,1,900,1200];
set(gcf, 'Position', sc);

for cc = 1:size(combinations_,1)*3

    subplot(9,3,cc)
    sgtitle({['Retrieval'];[]})

    contourf(MI.comod.parameters.PhaseFreqVector+MI.comod.parameters.PhaseFreq_BandWidth/2, ...
        MI.comod.parameters.AmpFreqVector+MI.comod.parameters.AmpFreq_BandWidth/2, squeeze(mean(toplot_{cc},3))',80,'lines','none')
    
    set(gca,'fontsize',6)
    
    x = xlim;
    y = ylim;    

%     xlim([x(1) 20])
%     ylim([y(1) 100])
%     xticks([2,5,10, 15, 20])
%     yticks([20, 40, 60, 80, 100])

    colormap parula
    b = colorbar;
    %clim([0.5*10^-3 2*10^-3])
    clims = clim();
    clim([clims(2)/10 clims(2)/2])
    set(b,'Ticks',[0.0001 clims(2)])


    set(gca,'fontsize',6)

    if cc<=3
        title(titles_{cc})
    end

    if cc == 1 || cc == 4 || cc == 7 || cc == 10 || cc == 13 || cc == 16 || cc == 19 || cc == 22 || cc == 25
        ylabel({[y_combination{ceil(cc/3)}];[];'Amplitude'})
    end

    if cc>=16
        xlabel('Phase (Hz)')
    end

end

clear ('toplot_','titles_','y_combination','sc','cc','x','y','b','clims')

%% Save

%newStr = regexprep(files.id.name,'.mat','_');
%newStr = files.id(ms).name(1:end-8);
newStr = id(1:end-8);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_MI.comodulogram');

% save figure
set(gcf,'renderer', 'painters');
exportgraphics(gcf,strcat(name,'_mean','.png'),'Resolution',300);

close all
clear('name','newStr','path') 

%% Graph comodulogram - all trials

figure
set(gcf,'color','w');
sc = [1,1,1800,600];
set(gcf, 'Position', sc);

comb = 7;

subplot(2,size(MI.comod.Comodulogram_CS{comb,1},3)+1,1)
sgtitle('habituation')

toplot_baseline = rescale(MI.comod.Comodulogram_baseline{comb,1});

contourf(MI.comod.parameters.PhaseFreqVector+MI.comod.parameters.PhaseFreq_BandWidth/2, MI.comod.parameters.AmpFreqVector+MI.comod.parameters.AmpFreq_BandWidth/2, toplot_baseline',80,'lines','none')
set(gca,'fontsize',12)
ylabel('Amplitude Frequency (Hz)')
xlabel('Phase Frequency (Hz)')
% xlim([3 12])
% ylim([35 140])

colormap parula
b = colorbar;
clims = clim();
clim([clims(2)/10 clims(2)/2])
set(b,'Ticks',[0.0001 clims(2)])

title ('Baseline')



for ii = 1:size(MI.comod.Comodulogram_CS{comb,1},3)

    subplot(2,CS_+1,ii+1)
    sgtitle('Retrieval')

    toplot = rescale(MI.comod.Comodulogram_CS{5,1}(:,:,ii));

    contourf(MI.comod.parameters.PhaseFreqVector+MI.comod.parameters.PhaseFreq_BandWidth/2, MI.comod.parameters.AmpFreqVector+MI.comod.parameters.AmpFreq_BandWidth/2, toplot',80,'lines','none')
    set(gca,'fontsize',12)
    ylabel('Amplitude Frequency (Hz)')
    xlabel('Phase Frequency (Hz)')
%     xlim([3 12])
%     ylim([35 140])
    
    colormap parula
    b = colorbar;
    clims = clim();
    clim([clims(2)/10 clims(2)/2])
    set(b,'Ticks',[0.0001 clims(2)])

    title ({['CS: ' num2str(ii)]})
end

for ii = 1:size(MI.comod.Comodulogram_ITI{comb,1},3)

    subplot(2,ITI_+1,ii+7)
    sgtitle('Retrieval')

    toplot = rescale(MI.comod.Comodulogram_ITI{5,1}(:,:,ii));

    contourf(MI.comod.parameters.PhaseFreqVector+MI.comod.parameters.PhaseFreq_BandWidth/2, MI.comod.parameters.AmpFreqVector+MI.comod.parameters.AmpFreq_BandWidth/2, toplot',80,'lines','none')
    set(gca,'fontsize',12)
    ylabel('Amplitude Frequency (Hz)')
    xlabel('Phase Frequency (Hz)')
%     xlim([3 12])
%     ylim([35 140])

    colormap parula
    b = colorbar;
    clims = clim();
    clim([clims(2)/10 clims(2)/2])
    set(b,'Ticks',[0.0001 clims(2)])

    title ({['ITI: ' num2str(ii)]})
end

clear('toplot_baseline','toplot') 

%% Save

%newStr = regexprep(files.id.name,'.mat','_');
%newStr = files.id(ms).name(1:end-8);
newStr = id(1:end-8);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_MI.comodulogram_trials');

% save data
save(name,'MI','-v7.3')

% save figure
set(gcf,'renderer', 'painters');
exportgraphics(gcf,strcat(name,'_trials', '.png'),'Resolution',300)

close all
clear('name','newStr','path') 

%% last update 04/04/2024 - 12:30
%  listening: Surrenderdorothy - IllLeavelUpToYou

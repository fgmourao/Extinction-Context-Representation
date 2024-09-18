% Final Plots Phase-amplitude Cross-frequency coupling measure

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  03/2024
% Last update: 03/2024

%% Load data retrieval

% Comodulogram

for ii = 1:size(MI.comod.parameters.ch_phase,1)

    data_2_plot_comod_baseline_retrieval{ii,1}(:,:,ms)   = MI.comod.Comodulogram_baseline{ii,1};
    data_2_plot_comod_CS_retrieval{ii,1}(:,:,:,ms)       = MI.comod.Comodulogram_CS{ii,1};
    data_2_plot_comod_ITI_retrieval{ii,1}(:,:,:,ms)      = MI.comod.Comodulogram_ITI{ii,1};

end


%% Averaging


if ms == length(files.FilesLoaded{1, 1})


    for ii = 1:size(MI.comod.parameters.ch_phase,1)

        data_2_plot_comod_baseline_retrieval{ii,1}(data_2_plot_comod_baseline_retrieval{ii,1} == 0) = NaN;
        data_2_plot_comod_CS_retrieval{ii,1}(data_2_plot_comod_CS_retrieval{ii,1} == 0) = NaN;
        data_2_plot_comod_ITI_retrieval{ii,1}(data_2_plot_comod_ITI_retrieval{ii,1} == 0) = NaN;

    end


    % Mean and SEM

    data_2_plot_comod_baseline_retrieval_Trial_mean = [];
    data_2_plot_comod_baseline_retrieval_Trial_std = [];
    data_2_plot_comod_CS_retrieval_Trial_mean = [];
    data_2_plot_comod_ITI_retrieval_Trial_mean = [];

    for ii = 1:size(MI.comod.parameters.ch_phase,1)

        data_2_plot_comod_baseline_retrieval_Trial_mean{ii,1}       = mean(data_2_plot_comod_baseline_retrieval{ii,1},3,'omitnan');
        data_2_plot_comod_baseline_retrieval_Trial_std{ii,1}        = std(data_2_plot_comod_baseline_retrieval{ii,1},[],3,'omitnan');

        data_2_plot_comod_CS_retrieval_Trial_mean{ii,1}             = max(mean(data_2_plot_comod_CS_retrieval{ii,1},4,'omitnan'),[],3,'omitnan');% - min(mean(data_2_plot_comod_CS_retrieval{ii,1}(:,:,[1 2],2:6),4,'omitnan'),[],3,'omitnan');
        data_2_plot_comod_ITI_retrieval_Trial_mean{ii,1}            = max(mean(data_2_plot_comod_ITI_retrieval{ii,1},4,'omitnan'),[],3,'omitnan');% - min(mean(data_2_plot_comod_ITI_retrieval{ii,1}(:,:,[1 2],2:6),4,'omitnan'),[],3,'omitnan');

%         data_2_plot_comod_CS_retrieval_Trial_mean{ii,1}             = rescale(squeeze(mean(data_2_plot_comod_CS_retrieval{ii,1}(:,:,1,:),4,'omitnan')));
%         data_2_plot_comod_ITI_retrieval_Trial_mean{ii,1}            = rescale(squeeze(mean(data_2_plot_comod_ITI_retrieval{ii,1}(:,:,1,:),4,'omitnan')));        

    end


    %% Graph comodulogram
    %  Events pars, avereged trials - Baseline, CS-Tones and ITI

    % All Possibles combinations
combinations_ = zeros(size(data.lfp{6,1},1)*3,2);
combinations_(1:2:6,:) = flip(nchoosek(1:size(data.lfp{6,1},1),2),2);
combinations_(2:2:6,:) = nchoosek(1:size(data.lfp{6,1},1),2);
combinations_(7:9,1) = 1:size(data.lfp{6,1},1);
combinations_(7:9,2) = 1:size(data.lfp{6,1},1);


    toplot_       = cell(max(max(combinations_)),size(combinations_,1));
    toplot_(1,:)  = data_2_plot_comod_baseline_retrieval_Trial_mean;
    toplot_(2,:)  = data_2_plot_comod_CS_retrieval_Trial_mean;
    toplot_(3,:)  = data_2_plot_comod_ITI_retrieval_Trial_mean;

%     for ii = 1:size(toplot_,2)
%         toplot_{2,ii} = (data_2_plot_comod_CS_retrieval_Trial_mean{ii,1} - data_2_plot_comod_baseline_retrieval_Trial_mean{ii,1})./data_2_plot_comod_baseline_retrieval_Trial_std{ii,1};
%     end

%%

% Frequency vector
PhaseFreqVector    = MI.comod.parameters.PhaseFreqVector;+MI.comod.parameters.PhaseFreq_BandWidth/2;
AmpFreqVector      = MI.comod.parameters.AmpFreqVector;+MI.comod.parameters.AmpFreq_BandWidth/2;


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

    for cc = 1:size(MI.comod.parameters.ch_phase,1)*3

        subplot(9,3,cc)
        sgtitle({['Retrieval'];[]})

        contourf(PhaseFreqVector, ...
            AmpFreqVector, toplot_{cc}',80,'lines','none')

        x = xlim;
        y = ylim;

        xlim([2.5 12])
        ylim([20 150])
        xticks([2,4,6,8,10,12])
        % ticks([20, 40, 60, 80, 100])

        colormap parula
        b = colorbar;

%         clims = clim();
%         clim([clims(2)/10 2*10^-3])

        %set(b,'Ticks',[0.03*10^-4 3*10^-4])


        set(gca,'fontsize',7)

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

    clear ('titles_','y_combination','sc','cc','x','y','b','clims')
    %% Save

    %newStr = regexprep(files.id.name,'.mat','_');
    %newStr = files.id(ms).name(1:end-8);
    newStr = id(1:end-8);

    path = '/Users/flavio/Desktop';
    %path = files.FilesLoaded{1, 1}.folder;

    name = strcat(path,'/',newStr,'_MI.comodulogram_retrieval');

    % save figure
    set(gcf,'renderer', 'painters');
    exportgraphics(gcf,strcat(name,'.png'),'Resolution',300)
    exportgraphics(gcf,strcat(name,'.eps'),'Resolution',300)
    
    %saveas(gcf,name,'png')

    close all
    clear('name','newStr','path')

end

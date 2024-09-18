% Final Plots Welch

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  01/2024
% Last update: 02/2024

% Extinction - CS-Trials without noise
% CSIT_extinction{1} = 1:5;      % MT1
% CSIT_extinction{2} = 1:4;      % MT3
% CSIT_extinction{3} = [1 2 4];  % MT4
% CSIT_extinction{4} = 1:5;      % MT5
% CSIT_extinction{5} = 1:5;      % MT6
% CSIT_extinction{6} = 1:5;      % MT7
% 
% final_plot_all_extinction(:,:,ms)  = mean(pw.stats_total_power{1,2}(:,:,CSIT_extinction{ms}),3);

% Retrieval - CS-Trials without noise
CSIT_retrieval{1} = 1:5;         % MT1
CSIT_retrieval{2} = [1 3 4 5];   % MT3
CSIT_retrieval{3} = [1];         % MT4
CSIT_retrieval{4} = 1:5;         % MT5
CSIT_retrieval{5} = 1:5;         % MT6
CSIT_retrieval{6} = 1:5;         % MT7

final_plot_all_retrieval(:,:,ms)  = mean(pw.stats_total_power{1,2}(:,:,CSIT_retrieval{ms}),3);


%% Total Power Normalization according to the Maren & Totty paper (pw.stats_total_power)


% Averaging all animals
%final_plot_extinction_mean  = smoothdata(mean(final_plot_all_extinction,3), 2 ,'gaussian',[6 2]);% Totty code. why?
final_plot_extinction_mean  = 10.*(mean(final_plot_all_extinction,3));
final_plot_extinction_SEM   = 10.*(std(final_plot_all_extinction,[],3)./size(final_plot_all_extinction,3));

%final_plot_retrieval_mean  = smoothdata(mean(final_plot_all_retrieval,3), 2 ,'gaussian',[6 2]); % Totty code. why?
final_plot_retrieval_mean  = 10.*(mean(final_plot_all_retrieval,3));
final_plot_retrieval_SEM   = 10.*(std(final_plot_all_retrieval,[],3)./size(final_plot_all_retrieval,3));

% Frequency range
steps         = diff(pw.full_trial.freq_CS_trials); % according to the fft time window
freq2plot_1   = 2:steps(1):12;
closestfreq_1 = dsearchn(pw.full_trial.freq_CS_trials,freq2plot_1');

rangef = pw.full_trial.freq_CS_trials(closestfreq_1);


%% Figure
figure
set(gcf,'color','w');
set(gcf, 'Position', get(0, 'Screensize'));

subplot(1,3,1)

% SEM shades
curve1 = final_plot_extinction_mean(1,:) + final_plot_extinction_SEM(1,:);
curve2 = final_plot_extinction_mean(1,:) - final_plot_extinction_SEM(1,:);
x2 = [rangef', fliplr(rangef')];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [.6, 0, 0],'FaceAlpha',0.2,'EdgeColor',[1 1 1]);

hold on

plot(rangef,final_plot_extinction_mean(1,:),'linew', 2,'Color',[.6, 0, 0])

% SEM shades
curve1 = final_plot_retrieval_mean(1,:) + final_plot_retrieval_SEM(1,:);
curve2 = final_plot_retrieval_mean(1,:) - final_plot_retrieval_SEM(1,:);
x2 = [rangef', fliplr(rangef')];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0, 0.7000, 0.9000],'FaceAlpha',0.2,'EdgeColor',[1 1 1]);

plot(rangef,final_plot_retrieval_mean(1,:),'linew', 2,'Color',[0, 0.7000, 0.9000])
xline(6)

xlabel('(Hz)','FontSize',9), ylabel('Normalized Power (A.U.)','FontSize',11)
xlim([2 12])

title ('PL')
ylim([0.001 0.15])
box off



subplot(1,3,2)

% SEM shades

curve1 = final_plot_extinction_mean(2,:) + final_plot_extinction_SEM(2,:);
curve2 = final_plot_extinction_mean(2,:) - final_plot_extinction_SEM(2,:);
x2 = [rangef', fliplr(rangef')];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [.6, 0, 0],'FaceAlpha',0.2,'EdgeColor',[1 1 1]);

hold on

plot(rangef,final_plot_extinction_mean(2,:),'linew', 2,'Color',[.6, 0, 0])

curve1 = final_plot_retrieval_mean(2,:) + final_plot_retrieval_SEM(2,:);
curve2 = final_plot_retrieval_mean(2,:) - final_plot_retrieval_SEM(2,:);
x2 = [rangef', fliplr(rangef')];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0, 0.7000, 0.9000],'FaceAlpha',0.2,'EdgeColor',[1 1 1]);

plot(rangef,final_plot_retrieval_mean(2,:),'linew', 2,'Color',[0, 0.7000, 0.9000])
xline(6)
xlabel('(Hz)','FontSize',9), ylabel('Normalized Power (A.U.)','FontSize',11)
xlim([2 12])
title ('IL')
ylim([0.001 0.15])
box off



subplot(1,3,3)

curve1 = final_plot_extinction_mean(3,:) + final_plot_extinction_SEM(3,:);
curve2 = final_plot_extinction_mean(3,:) - final_plot_extinction_SEM(3,:);
x2 = [rangef', fliplr(rangef')];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [.6, 0, 0],'FaceAlpha',0.2,'EdgeColor',[1 1 1]);

hold on

plot(rangef,final_plot_extinction_mean(3,:),'linew', 2,'Color',[.6, 0, 0])

curve1 = final_plot_retrieval_mean(3,:) + final_plot_retrieval_SEM(3,:);
curve2 = final_plot_retrieval_mean(3,:) - final_plot_retrieval_SEM(3,:);
x2 = [rangef', fliplr(rangef')];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0, 0.7000, 0.9000],'FaceAlpha',0.2,'EdgeColor',[1 1 1]);

plot(rangef,final_plot_retrieval_mean(3,:),'linew', 2,'Color',[0, 0.7000, 0.9000])
xline(6)
xlabel('(Hz)','FontSize',9), ylabel('Normalized Power (A.U.)','FontSize',11)
xlim([2 12])
title ('HPC')
ylim([0.001 0.15])
box off







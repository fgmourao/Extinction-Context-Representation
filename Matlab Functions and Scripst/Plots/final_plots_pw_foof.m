% Final Plots Welch

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  01/2024
% Last update: 02/2024

%final_plot_all_extinction(:,:,ms)  = pw.over_f.stats_total_power_mean{1, 2};
% final_plot_all_retrieval(:,:,ms)   = pw.over_f.stats_total_power_mean{1, 2};


%% Total Power Normalization according to the Maren & Totty paper (pw.stats_total_power)


% Averaging all animals
%final_plot_extinction_mean  = smoothdata(mean(final_plot_all_extinction,3), 2 ,'gaussian',[6 2]);% Totty code. why?
final_plot_extinction_mean  = mean(final_plot_all_extinction,3);
final_plot_extinction_SEM   = std(final_plot_all_extinction,[],3)./size(final_plot_all_extinction,3);

%final_plot_retrieval_mean  = smoothdata(mean(final_plot_all_retrieval,3), 2 ,'gaussian',[6 2]); % Totty code. why?
final_plot_retrieval_mean  = mean(final_plot_all_retrieval,3);
final_plot_retrieval_SEM   = std(final_plot_all_retrieval,[],3)./size(final_plot_all_retrieval,3);

% Frequency range
theta_range = [2 12];
closestfreq_1  = find(theta_range(1)<pw.over_f.CS_Trials(1, 1).Analise.freqs & pw.over_f.CS_Trials(1, 1).Analise.freqs<theta_range(2));
rangef = pw.over_f.CS_Trials(1, 1).Analise.freqs(closestfreq_1);


% Figure
figure
set(gcf,'color','w');
set(gcf, 'Position', get(0, 'Screensize'));

subplot(1,3,1)

% SEM shades
curve1 = final_plot_extinction_mean(1,:) + final_plot_extinction_SEM(1,:);
curve2 = final_plot_extinction_mean(1,:) - final_plot_extinction_SEM(1,:);
x2 = [rangef, fliplr(rangef)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [.6, 0, 0],'FaceAlpha',0.2,'EdgeColor',[1 1 1]);

hold on

plot(rangef,final_plot_extinction_mean(1,:),'linew', 2,'Color',[.6, 0, 0])

% SEM shades
curve1 = final_plot_retrieval_mean(1,:) + final_plot_retrieval_SEM(1,:);
curve2 = final_plot_retrieval_mean(1,:) - final_plot_retrieval_SEM(1,:);
x2 = [rangef, fliplr(rangef)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0, 0.7000, 0.9000],'FaceAlpha',0.2,'EdgeColor',[1 1 1]);

plot(rangef,final_plot_retrieval_mean(1,:),'linew', 2,'Color',[0, 0.7000, 0.9000])
xline(6)

xlabel('(Hz)','FontSize',9), ylabel('Normalized Power (A.U.)','FontSize',11)
xlim([2 12])

title ('PL')
ylim([0.001 0.10])
box off



subplot(1,3,2)

% SEM shades

curve1 = final_plot_extinction_mean(2,:) + final_plot_extinction_SEM(2,:);
curve2 = final_plot_extinction_mean(2,:) - final_plot_extinction_SEM(2,:);
x2 = [rangef, fliplr(rangef)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [.6, 0, 0],'FaceAlpha',0.2,'EdgeColor',[1 1 1]);

hold on

plot(rangef,final_plot_extinction_mean(2,:),'linew', 2,'Color',[.6, 0, 0])

curve1 = final_plot_retrieval_mean(2,:) + final_plot_retrieval_SEM(2,:);
curve2 = final_plot_retrieval_mean(2,:) - final_plot_retrieval_SEM(2,:);
x2 = [rangef, fliplr(rangef)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0, 0.7000, 0.9000],'FaceAlpha',0.2,'EdgeColor',[1 1 1]);

plot(rangef,final_plot_retrieval_mean(2,:),'linew', 2,'Color',[0, 0.7000, 0.9000])
xline(6)
xlabel('(Hz)','FontSize',9), ylabel('Normalized Power (A.U.)','FontSize',11)
xlim([2 12])
title ('IL')
ylim([0.001 0.10])
box off



subplot(1,3,3)

curve1 = final_plot_extinction_mean(3,:) + final_plot_extinction_SEM(3,:);
curve2 = final_plot_extinction_mean(3,:) - final_plot_extinction_SEM(3,:);
x2 = [rangef, fliplr(rangef)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [.6, 0, 0],'FaceAlpha',0.2,'EdgeColor',[1 1 1]);

hold on

plot(rangef,final_plot_extinction_mean(3,:),'linew', 2,'Color',[.6, 0, 0])

curve1 = final_plot_retrieval_mean(3,:) + final_plot_retrieval_SEM(3,:);
curve2 = final_plot_retrieval_mean(3,:) - final_plot_retrieval_SEM(3,:);
x2 = [rangef, fliplr(rangef)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0, 0.7000, 0.9000],'FaceAlpha',0.2,'EdgeColor',[1 1 1]);

plot(rangef,final_plot_retrieval_mean(3,:),'linew', 2,'Color',[0, 0.7000, 0.9000])
xline(6)
xlabel('(Hz)','FontSize',9), ylabel('Normalized Power (A.U.)','FontSize',11)
xlim([2 12])
title ('HPC')
ylim([0.001 0.10])
box off







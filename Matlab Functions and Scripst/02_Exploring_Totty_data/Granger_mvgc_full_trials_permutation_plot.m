%% Extinction

steps                           = diff(mvgc.parameters.freqs); % according to the fft time window
mvgc.parameters.frex_2_12Hz     = 2:steps(1):12;
mvgc.parameters.frex_idx_2_12Hz = dsearchn(mvgc.parameters.freqs,mvgc.parameters.frex_2_12Hz');
freq_v = mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz);q


Granger_mvgc_permutation_baseline_mean = mean(cat(5,Granger_mvgc_permutation_baseline{1,:}),5);
Granger_mvgc_permutation_baseline_SEM  = std(cat(5,Granger_mvgc_permutation_baseline{1,:}),[],5)./size(Granger_mvgc_permutation_baseline,2);

Granger_mvgc_permutation_first_CS_Tone_mean = mean(cat(5,Granger_mvgc_permutation_first_CS_Tone{1,:}),5);
Granger_mvgc_permutation_first_CS_Tone_SEM  = std(cat(5,Granger_mvgc_permutation_first_CS_Tone{1,:}),[],5)./size(Granger_mvgc_permutation_first_CS_Tone,2);

Granger_mvgc_permutation_last_CS_Tone_mean = mean(cat(5,Granger_mvgc_permutation_last_CS_Tone{1,:}),5);
Granger_mvgc_permutation_last_CS_Tone_SEM  = std(cat(5,Granger_mvgc_permutation_last_CS_Tone{1,:}),[],5)./size(Granger_mvgc_permutation_last_CS_Tone,2);

%% Retrieval

steps                           = diff(mvgc.parameters.freqs); % according to the fft time window
mvgc.parameters.frex_2_12Hz     = 2:steps(1):12;
mvgc.parameters.frex_idx_2_12Hz = dsearchn(mvgc.parameters.freqs,mvgc.parameters.frex_2_12Hz');
freq_v = mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz);


Granger_mvgc_permutation_baseline_mean = mean(cat(5,Granger_mvgc_permutation_baseline{1,:}),5);
Granger_mvgc_permutation_baseline_SEM  = std(cat(5,Granger_mvgc_permutation_baseline{1,:}),[],5)./size(Granger_mvgc_permutation_baseline,2);

Granger_mvgc_permutation_CS_Tone_mean = mean(cat(5,Granger_mvgc_permutation_CS_Tone{1,:}),5);
Granger_mvgc_permutation_CS_Tone_SEM  = std(cat(5,Granger_mvgc_permutation_CS_Tone{1,:}),[],5)./size(Granger_mvgc_permutation_CS_Tone,2);


%%

figure
set(gcf,'color','w');
sc = [1,1,960,1200];
set(gcf, 'Position', sc);
sgtitle({'Pairwise-conditional Granger causality - frequency domain.';''});

subplot(2,3,1)
boundedline(freq_v',squeeze(mean(Granger_mvgc_permutation_baseline_mean(:,1,3,mvgc.parameters.frex_idx_2_12Hz),1)),squeeze(mean(Granger_mvgc_permutation_baseline_SEM(:,1,3,mvgc.parameters.frex_idx_2_12Hz),1)),'LineWidth',2,'color',[.6 .6 .6],'transparency',.4);
hold on
boundedline(freq_v',squeeze(mean(Granger_mvgc_permutation_baseline_mean(:,3,1,mvgc.parameters.frex_idx_2_12Hz),1)),squeeze(mean(Granger_mvgc_permutation_baseline_SEM(:,3,1,mvgc.parameters.frex_idx_2_12Hz),1)),'LineWidth',2,'color',[.3 .3 .3],'transparency',.4);
ylabel({'Granger causality'},'FontWeight','bold','FontSize',12)
xlabel('(Hz)','FontSize',12)
title("baseline")
legend('IL --> HPC','','HPC --> IL','')
legend('boxoff')
legend('FontSize',8)
ylim([0 0.2])

subplot(2,3,2)
boundedline(freq_v',squeeze(mean(Granger_mvgc_permutation_CS_Tone_mean(:,1,3,mvgc.parameters.frex_idx_2_12Hz),1)),squeeze(mean(Granger_mvgc_permutation_first_CS_Tone_SEM(:,1,3,mvgc.parameters.frex_idx_2_12Hz),1)),'LineWidth',2,'color',[.6 .6 .6],'transparency',.4);
hold on
boundedline(freq_v',squeeze(mean(Granger_mvgc_permutation_CS_Tone_mean(:,3,1,mvgc.parameters.frex_idx_2_12Hz),1)),squeeze(mean(Granger_mvgc_permutation_first_CS_Tone_SEM(:,3,1,mvgc.parameters.frex_idx_2_12Hz),1)),'LineWidth',2,'color',[.3 .3 .3],'transparency',.4);
ylabel({'Granger causality'},'FontWeight','bold','FontSize',12)
xlabel('(Hz)','FontSize',12)
title("First 10 Trials")
legend('IL --> HPC','','HPC --> IL','')
legend('boxoff')
legend('FontSize',8)
ylim([0 0.2])


subplot(2,3,4)
boundedline(freq_v',squeeze(mean(Granger_mvgc_permutation_baseline_mean(:,2,3,mvgc.parameters.frex_idx_2_12Hz),1)),squeeze(mean(Granger_mvgc_permutation_baseline_SEM(:,2,3,mvgc.parameters.frex_idx_2_12Hz),1)),'LineWidth',2,'color',[.6 .6 .6],'transparency',.4);
hold on
boundedline(freq_v',squeeze(mean(Granger_mvgc_permutation_baseline_mean(:,3,2,mvgc.parameters.frex_idx_2_12Hz),1)),squeeze(mean(Granger_mvgc_permutation_baseline_SEM(:,3,2,mvgc.parameters.frex_idx_2_12Hz),1)),'LineWidth',2,'color',[.3 .3 .3],'transparency',.4);
ylabel({'Granger causality'},'FontWeight','bold','FontSize',12)
xlabel('(Hz)','FontSize',12)
title("baseline")
legend('IL --> HPC','','HPC --> IL','')
legend('boxoff')
legend('FontSize',8)
ylim([0 0.2])

subplot(2,3,5)
boundedline(freq_v',squeeze(mean(Granger_mvgc_permutation_CS_Tone_mean(:,2,3,mvgc.parameters.frex_idx_2_12Hz),1)),squeeze(mean(Granger_mvgc_permutation_first_CS_Tone_SEM(:,2,3,mvgc.parameters.frex_idx_2_12Hz),1)),'LineWidth',2,'color',[.6 .6 .6],'transparency',.4);
hold on
boundedline(freq_v',squeeze(mean(Granger_mvgc_permutation_CS_Tone_mean(:,3,2,mvgc.parameters.frex_idx_2_12Hz),1)),squeeze(mean(Granger_mvgc_permutation_first_CS_Tone_SEM(:,3,2,mvgc.parameters.frex_idx_2_12Hz),1)),'LineWidth',2,'color',[.3 .3 .3],'transparency',.4);
ylabel({'Granger causality'},'FontWeight','bold','FontSize',12)
xlabel('(Hz)','FontSize',12)
title("First 10 Trials")
legend('IL --> HPC','','HPC --> IL','')
legend('boxoff')
legend('FontSize',8)
ylim([0 0.2])



%%
%zscore(squeeze(mean(Granger_mvgc_(:,1,3,:),4)));

zscore(squeeze(mean(Granger_mvgc_(:,3,1,:),4)),[],1);
histogram(ans,60)
histfit(ans)
hold on

z_observed = (squeeze(mean(mvgc.F_spect{1,1}(3,1,mvgc.parameters.frex_idx_6_8Hz),3))...
    - mean(mean(Granger_mvgc_(:,3,1,:),4),1))/...
    std(mean(Granger_mvgc_(:,3,1,:),4),[],1);% z real(observed) value


p_value = 0.5 * erfc(z_observed  ./ sqrt(2)); % p value. Similar to Matlab function: normcdf(-z) two-tailed

stem(z_observed,50,'linew',2);
%xlim([-3.5 6])

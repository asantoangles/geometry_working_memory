%% GERE project
% replay

% path to data
if isfolder('/path_to_local')
    path_data = '/path_to_local/results/preprocessing/M3';
    path_inputs = ['/path_to_local/results/mvpa/' folder_decoding];
    path_results = ['/path_to_local/results/replay/' folder_replay];
    path_sequences = '/path_to_local/results/replay/sequences';
    classifiers = {'l2_e-2'};
    inputs = {'bins'};
    state_spaces = {'prob'};
    addpath('/path_to_local/scripts/replay/utilities')
end

number_stimuli = 8;
folders = {'task_stim2delay_balanced_trial_predictions_averaged'...
           'localizer2delay_balanced_trial_predictions_averaged'};
folder_output = 'TDLM_sequence';
time_windows = {'whole_delay'};

%% replay

for folder_i = 1:length(folders)

    predictions_folder = folders{folder_i};

    for time_window_i = 1:length(time_windows)
        
        for state_i = 1:length(state_spaces)

            for input_i = 1:length(inputs)
    
                if strcmp(inputs{input_i}, 'timepoints')
                    nlags = 500; % number of different lags (500 medians that shifts include [1...500] timepoints of the period of interest (delay period)
                elseif strcmp(inputs{input_i}, 'bins')
                    nlags = 100;
                end
                
                if strcmp(inputs{input_i}, 'timepoints')
                    if strcmp(time_windows{time_window_i}, 'whole_delay')
                        time_window = 1:4000;
                    elseif strcmp(time_windows{time_window_i}, 'early_delay')
                        time_window = 1:2000;
                    elseif strcmp(time_windows{time_window_i}, 'late_delay')
                        time_window = 2001:4000;
                    end
                elseif strcmp(inputs{input_i}, 'bins')
                    if strcmp(time_windows{time_window_i}, 'whole_delay')
                        time_window = 1:791;
                    elseif strcmp(time_windows{time_window_i}, 'early_delay')
                        time_window = 1:400;
                    elseif strcmp(time_windows{time_window_i}, 'late_delay')
                        time_window = 401:791;
                    end
                end
    
                for class_i = 1:length(classifiers)
        
                    classifier = classifiers{class_i};

                    %% pool data across subjects

                    for sub_i = 1:length(subjects)
                    
                        subject = subjects(sub_i);
                                        
                        % set paths
                        if subject < 10
                            subject_ID = ['sub_0' num2str(subject)];
                            subjectID = ['sub0' num2str(subject)];
                        else
                            subject_ID = ['sub_' num2str(subject)];
                            subjectID = ['sub' num2str(subject)];
                        end
                            
                        load([path_results '/' subject_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/Msf_l3.mat']);
                        load([path_results '/' subject_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/Msf_l4.mat']);
                        load([path_results '/' subject_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/Msb_l3.mat']);
                        load([path_results '/' subject_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/Msb_l4.mat']);                            
                        load([path_results '/' subject_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/Msf_l34.mat']);
                        load([path_results '/' subject_ID '/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/Msb_l34.mat']);   

                        Msf_l3 = Msf_l3(:,:,1:(end/2));
                        Msb_l3 = Msb_l3(:,:,1:(end/2));
                        Msf_l4 = Msf_l4(:,:,1:(end/2));
                        Msb_l4 = Msb_l4(:,:,1:(end/2));

                        if sub_i == 1

                            Msf_l3_allsubjects = Msf_l3;
                            Msb_l3_allsubjects = Msb_l3;
                            Msf_l4_allsubjects = Msf_l4;
                            Msb_l4_allsubjects = Msb_l4;
                            Msf_l34_allsubjects = Msf_l34;
                            Msb_l34_allsubjects = Msb_l34;

                        else

                            Msf_l3_allsubjects(end+1,:,:) = Msf_l3;
                            Msb_l3_allsubjects(end+1,:,:) = Msb_l3;
                            Msf_l4_allsubjects(end+1,:,:) = Msf_l4;
                            Msb_l4_allsubjects(end+1,:,:) = Msb_l4;
                            Msf_l34_allsubjects(end+1,:,:) = Msf_l34;
                            Msb_l34_allsubjects(end+1,:,:) = Msb_l34;

                        end

                    end

                    %% pool individual figures
                    
                    if ~isfile([path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                        folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures/sequenceness_l4_500.jpeg'])

                        %%% length 3
                        figure
    
                        for sub_ses_i = 1:size(Msf_l3_allsubjects)
                                            
                            subplot(4,5,sub_ses_i);
                            plot(squeeze(Msf_l3_allsubjects(sub_ses_i,1,:)), 'LineWidth', 1); hold on
                            plot(squeeze(Msb_l3_allsubjects(sub_ses_i,1,:)), 'LineWidth', 1)    
                            yline(max(max(Msf_l3_allsubjects(sub_ses_i,2:end,:))),'--b', 'LineWidth', 0.5);
                            yline(min(min(Msf_l3_allsubjects(sub_ses_i,2:end,:))),'--b', 'LineWidth', 0.5);
                            yline(max(max(Msb_l3_allsubjects(sub_ses_i,2:end,:))),'--r', 'LineWidth', 0.5);
                            yline(min(min(Msb_l3_allsubjects(sub_ses_i,2:end,:))),'--r', 'LineWidth', 0.5);
                            yline(0,'-k', 'LineWidth', 0.5);
    
                            set(gca, 'XTick', [], 'YTick', []);
                            set(gcf, 'Color', 'white');
                            set(gca, 'Color', 'white', 'Box', 'on');
    
                        end
    
                        saveas(gcf, [path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                            folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures/sequenceness_l3_500.fig']);
                        fig = openfig([path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                            folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures/sequenceness_l3_500.fig']);
                        saveas(fig, [path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                            folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures/sequenceness_l3_500.jpeg']);
    
                        system(['rm ' path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                            folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures/sequenceness_l3_500.fig']);
    
                        close all
    
                        %%% length 4
                        figure
    
                        for sub_ses_i = 1:size(Msf_l3_allsubjects)
                                            
                            subplot(4,5,sub_ses_i);
                            plot(squeeze(Msf_l4_allsubjects(sub_ses_i,1,:)), 'LineWidth', 1); hold on
                            plot(squeeze(Msb_l4_allsubjects(sub_ses_i,1,:)), 'LineWidth', 1)    
                            yline(max(max(Msf_l4_allsubjects(sub_ses_i,2:end,:))),'--b', 'LineWidth', 0.5);
                            yline(min(min(Msf_l4_allsubjects(sub_ses_i,2:end,:))),'--b', 'LineWidth', 0.5);
                            yline(max(max(Msb_l4_allsubjects(sub_ses_i,2:end,:))),'--r', 'LineWidth', 0.5);
                            yline(min(min(Msb_l4_allsubjects(sub_ses_i,2:end,:))),'--r', 'LineWidth', 0.5);
                            yline(0,'-k', 'LineWidth', 0.5);
    
                            set(gca, 'XTick', [], 'YTick', []);
                            set(gcf, 'Color', 'white');
                            set(gca, 'Color', 'white', 'Box', 'on');
    
                        end
    
                        saveas(gcf, [path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                            folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures/sequenceness_l4_500.fig']);
                        fig = openfig([path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                            folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures/sequenceness_l4_500.fig']);
                        saveas(fig, [path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                            folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures/sequenceness_l4_500.jpeg']);
    
                        system(['rm ' path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                            folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures/sequenceness_l4_500.fig']);
    
                        close all
    
                        %%% length all
                        figure
    
                        for sub_ses_i = 1:size(Msf_l34_allsubjects)
                                            
                            subplot(4,5,sub_ses_i);
                            plot(squeeze(Msf_l34_allsubjects(sub_ses_i,1,:)), 'LineWidth', 1); hold on
                            plot(squeeze(Msb_l34_allsubjects(sub_ses_i,1,:)), 'LineWidth', 1)    
                            yline(max(max(Msf_l34_allsubjects(sub_ses_i,2:end,:))),'--b', 'LineWidth', 0.5);
                            yline(min(min(Msf_l34_allsubjects(sub_ses_i,2:end,:))),'--b', 'LineWidth', 0.5);
                            yline(max(max(Msb_l34_allsubjects(sub_ses_i,2:end,:))),'--r', 'LineWidth', 0.5);
                            yline(min(min(Msb_l34_allsubjects(sub_ses_i,2:end,:))),'--r', 'LineWidth', 0.5);
                            yline(0,'-k', 'LineWidth', 0.5);
    
                            set(gca, 'XTick', [], 'YTick', []);
                            set(gcf, 'Color', 'white');
                            set(gca, 'Color', 'white', 'Box', 'on');
    
                        end
    
                        saveas(gcf, [path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                            folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures/sequenceness_l34.fig']);
                        fig = openfig([path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                            folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures/sequenceness_l34.fig']);
                        saveas(fig, [path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                            folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures/sequenceness_l34.jpeg']);
    
                        system(['rm ' path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                            folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures/sequenceness_l34.fig']);
    
                        close all

                    end

                    %% figure - average across subjects

                    if ~isfile([path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                        folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures/sequenceness_500.jpeg'])

                        figure

                        % plot l3
                        subplot(3,1,1);
                        plot(squeeze(median(Msf_l3_allsubjects(:,1,:),1)), 'LineWidth', 2); hold on
                        plot(squeeze(median(Msb_l3_allsubjects(:,1,:),1)), 'LineWidth', 2)    
                        yline(median(max(Msf_l3_allsubjects(:,2:end,:),[],2:3)),'--b', 'LineWidth', 1.5);
                        yline(median(min(Msf_l3_allsubjects(:,2:end,:),[],2:3)),'--b', 'LineWidth', 1.5);
                        yline(median(max(Msb_l3_allsubjects(:,2:end,:),[],2:3)),'--r', 'LineWidth', 1.5);
                        yline(median(min(Msb_l3_allsubjects(:,2:end,:),[],2:3)),'--r', 'LineWidth', 1.5);
                        title('length 3');

                        % plot l4
                        subplot(3,1,2);
                        plot(squeeze(median(Msf_l4_allsubjects(:,1,:),1)), 'LineWidth', 2); hold on
                        plot(squeeze(median(Msb_l4_allsubjects(:,1,:),1)), 'LineWidth', 2)    
                        yline(median(max(Msf_l4_allsubjects(:,2:end,:),[],2:3)),'--b', 'LineWidth', 1.5);
                        yline(median(min(Msf_l4_allsubjects(:,2:end,:),[],2:3)),'--b', 'LineWidth', 1.5);
                        yline(median(max(Msb_l4_allsubjects(:,2:end,:),[],2:3)),'--r', 'LineWidth', 1.5);
                        yline(median(min(Msb_l4_allsubjects(:,2:end,:),[],2:3)),'--r', 'LineWidth', 1.5);
                        title('length 4');

                        % plot all lengths
                        subplot(3,1,3);
                        plot(squeeze(median(Msf_l34_allsubjects(:,1,:),1)), 'LineWidth', 2); hold on
                        plot(squeeze(median(Msb_l34_allsubjects(:,1,:),1)), 'LineWidth', 2)    
                        yline(median(max(Msf_l34_allsubjects(:,2:end,:),[],2:3)),'--b', 'LineWidth', 1.5);
                        yline(median(min(Msf_l34_allsubjects(:,2:end,:),[],2:3)),'--b', 'LineWidth', 1.5);
                        yline(median(max(Msb_l34_allsubjects(:,2:end,:),[],2:3)),'--r', 'LineWidth', 1.5);
                        yline(median(min(Msb_l34_allsubjects(:,2:end,:),[],2:3)),'--r', 'LineWidth', 1.5);
                        title('all lengths');

                        lgd = legend('forward', 'backward', 'Location', 'southoutside');
                        lgd.Position = [0.5 0.02 0.05 0.05];
                        lgd.Position(2) = lgd.Position(2) - 0.01;

                        saveas(gcf, [path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                            folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures/sequenceness_500.fig']);
                        fig = openfig([path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                            folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures/sequenceness_500.fig']);
                        saveas(fig, [path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                            folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures/sequenceness_500.jpeg']);

                        system(['rm ' path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                            folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i} '/figures/sequenceness_500.fig']);
                                 
                    end

                    %% figure single - average across subjects
    
                    % plot l3

                    null_Msf_thr_pos = median(max(Msf_l3_allsubjects(:,2:end,:),[],2:3));
                    null_Msf_thr_neg = median(min(Msf_l3_allsubjects(:,2:end,:),[],2:3));
                    null_Msb_thr_pos = median(max(Msb_l3_allsubjects(:,2:end,:),[],2:3));
                    null_Msb_thr_neg = median(min(Msb_l3_allsubjects(:,2:end,:),[],2:3));

                    % avoid full overlap of threshold lines
                    if round(null_Msf_thr_pos,3) == round(null_Msb_thr_pos,3)
                        null_Msf_thr_pos = null_Msf_thr_pos + 0.001;
                    end

                    if round(null_Msf_thr_neg,3) == round(null_Msb_thr_neg,3)
                        null_Msf_thr_neg = null_Msf_thr_neg + 0.001;
                    end

                    figure
                    plot((1:nlags)*5, squeeze(median(Msf_l3_allsubjects(:,1,:),1)), 'Color', [0 0 0.8], 'LineWidth', 4); hold on
                    plot((1:nlags)*5, squeeze(median(Msb_l3_allsubjects(:,1,:),1)), 'Color', [0.8 0 0], 'LineWidth', 4)    
                    yline(null_Msf_thr_pos,'--', 'Color', [0 0 0.8], 'LineWidth', 3);
                    yline(null_Msf_thr_neg,'--', 'Color', [0 0 0.8], 'LineWidth', 3);
                    yline(null_Msb_thr_pos,'--', 'Color', [0.8 0 0], 'LineWidth', 3);
                    yline(null_Msb_thr_neg,'--', 'Color', [0.8 0 0], 'LineWidth', 3);
    %                     title('length 3');
    %                     ylabel('sequenceness');
    %                     xlabel('time lag (ms)');
                    yticks(-0.06:0.02:0.06);
                    xticks([]);

                    if strcmp(predictions_folder, 'localizer2delay_balanced_trial_predictions_averaged')
                        ylabel('sequenceness');
                        xlabel('time lag (ms)');
                        xticks(0:100:500);
                    end
    
                    % Increase the font size of the tick labels
                    ax = gca; % Get current axes                            
                    ax.FontSize = 20;
                    
                    saveas(gcf, [path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                        folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i}  '/figures/sequenceness_length3.fig']);
                    fig = openfig([path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                        folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i}  '/figures/sequenceness_length3.fig']);
                    saveas(fig, [path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                        folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i}  '/figures/sequenceness_length3.jpeg']);
    
                    system(['rm ' path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                        folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i}  '/figures/sequenceness_length3.fig']);
    
    
                    % plot l4

                    null_Msf_thr_pos = median(max(Msf_l4_allsubjects(:,idx_null_l4,:),[],2:3));
                    null_Msf_thr_neg = median(min(Msf_l4_allsubjects(:,idx_null_l4,:),[],2:3));
                    null_Msb_thr_pos = median(max(Msb_l4_allsubjects(:,idx_null_l4,:),[],2:3));
                    null_Msb_thr_neg = median(min(Msb_l4_allsubjects(:,idx_null_l4,:),[],2:3));

                    % avoid full overlap of threshold lines
                    if round(null_Msf_thr_pos,3) == round(null_Msb_thr_pos,3)
                        null_Msf_thr_pos = null_Msf_thr_pos + 0.001;
                    end

                    if round(null_Msf_thr_neg,3) == round(null_Msb_thr_neg,3)
                        null_Msf_thr_neg = null_Msf_thr_neg + 0.001;
                    end

                    figure
                    plot((1:nlags)*5, squeeze(median(Msf_l4_allsubjects(:,1,:),1)), 'Color', [0 0 0.8], 'LineWidth', 4); hold on
                    plot((1:nlags)*5, squeeze(median(Msb_l4_allsubjects(:,1,:),1)), 'Color', [0.8 0 0], 'LineWidth', 4)    
                    yline(null_Msf_thr_pos,'--', 'Color', [0 0 0.8], 'LineWidth', 3);
                    yline(null_Msf_thr_neg,'--', 'Color', [0 0 0.8], 'LineWidth', 3);
                    yline(null_Msb_thr_pos,'--', 'Color', [0.8 0 0], 'LineWidth', 3);
                    yline(null_Msb_thr_neg,'--', 'Color', [0.8 0 0], 'LineWidth', 3);
    %                     title('length 3');
    %                     ylabel('sequenceness');
    %                     xlabel('time lag (ms)');
                    yticks(-0.06:0.02:0.06);
                    xticks([]);
    
                    % Increase the font size of the tick labels
                    ax = gca; % Get current axes                            
                    ax.FontSize = 20;
                    
                    saveas(gcf, [path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                        folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i}  '/figures/sequenceness_length4.fig']);
                    fig = openfig([path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                        folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i}  '/figures/sequenceness_length4.fig']);
                    saveas(fig, [path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                        folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i}  '/figures/sequenceness_length4.jpeg']);
    
                    system(['rm ' path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                        folder_output '/' predictions_folder '/' state_spaces{state_i} '/' time_windows{time_window_i}  '/figures/sequenceness_length4.fig']);

                end

            end
        
        end
    
    end
    
end

disp('Analysis done');
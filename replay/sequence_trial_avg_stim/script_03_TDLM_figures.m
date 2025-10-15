%% GERE project
% replay

% path to data
if isfolder('/path_to_local')
    path_data = '/path_to_local/results/preprocessing/M3';
    path_inputs = ['/path_to_local/results/mvpa/' folder_decoding];
    path_results = ['/path_to_local/results/replay/' folder_replay];
    path_sequences = '/path_to_local/results/replay/sequences';
    addpath('/path_to_local/scripts/replay/utilities')
end

classifiers = {'l2_e-1'};
inputs = {'bins'};
state_spaces = {'prob'};

number_stimuli = 8;
folders = {'localizer2task_stim_trial_predictions_averaged' 'task_stim_balanced2task_stim_trial_predictions_averaged'};
folder_output = 'TDLM_sequence';

%% stats - group-level

for folder_i = 1:length(folders)

    predictions_folder = folders{folder_i};

    disp(predictions_folder)

    try
    
        for state_i = 1:length(state_spaces)
        
            %% set input   
            for input_i = 1:length(inputs)
        
                if strcmp(inputs{input_i}, 'timepoints')
                    nlags = 500; % this variable can be shorter than 500, but no longer
                elseif strcmp(inputs{input_i}, 'bins')
                    nlags = 100;
                end
        
                %% loop across classifiers
                for class_i = 1:length(classifiers)
        
                    classifier = classifiers{class_i};
        
                    if isfile([path_results '/group_results/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/figures/data_individual_plots.mat'])
                
                        load([path_results '/group_results/' inputs{input_i} '/' classifier '/' folder_output '/' predictions_folder '/' state_spaces{state_i} '/figures/data_individual_plots.mat']);

                        if ~isfile([path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                                folder_output '/' predictions_folder '/' state_spaces{state_i} '/figures/Msfb_length4.jpeg'])
        
                            %% all lengths
                
                            % average across subjects/sessions
                            null_Msz_thr_pos = median(null_Msz_thr_pos_allsub);
                            null_Msz_thr_neg = median(null_Msz_thr_neg_allsub);
                            null_Msf_thr_pos = median(null_Msf_thr_pos_allsub);
                            null_Msf_thr_neg = median(null_Msf_thr_neg_allsub);
                            null_Msb_thr_pos = median(null_Msb_thr_pos_allsub);
                            null_Msb_thr_neg = median(null_Msb_thr_neg_allsub);

                            % avoid full overlap of threshold lines
                            if round(null_Msf_thr_pos,3) == round(null_Msb_thr_pos,3)
                                null_Msf_thr_pos = null_Msf_thr_pos + 0.005;
                            end

                            if round(null_Msf_thr_neg,3) == round(null_Msb_thr_neg,3)
                                null_Msf_thr_neg = null_Msf_thr_neg + 0.005;
                            end
            
                            
                            Msz = median(Msz_allsub, 2);
                            Msf = median(Msf_allsub, 2);
                            Msb = median(Msb_allsub, 2);

                            % plot Msf and Msb
                            figure;
                            plot((1:nlags)*5, Msf, 'Color', [0 0 0.8], 'LineWidth', 4); hold on;
                            plot((1:nlags)*5, Msb, 'Color', [0.8 0 0], 'LineWidth', 4); hold on;

                            yline(null_Msf_thr_pos,'--', 'Color', [0 0 0.8], 'LineWidth', 3);
                            yline(null_Msf_thr_neg,'--', 'Color', [0 0 0.8], 'LineWidth', 3);
                            yline(null_Msb_thr_pos,'--', 'Color', [0.8 0 0], 'LineWidth', 3);
                            yline(null_Msb_thr_neg,'--', 'Color', [0.8 0 0], 'LineWidth', 3);
%                             title(['Forward and Backward all lengths ' inputs{input_i} ' ' classifier ' ' state_spaces{state_i}])
%                             legend('forward', 'backward', 'Location', 'northwest');
%                             ylabel('sequenceness');
%                             xlabel('time lag (ms)');
                            yticks([-0.1 0 0.1 0.2]);
                            xticks([]);

                            % Increase the font size of the tick labels
                            ax = gca; % Get current axes                            
                            ax.FontSize = 20;

                
                            saveas(gcf, [path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                                folder_output '/' predictions_folder '/' state_spaces{state_i} '/figures/Msfb_lengthall.fig']);
                            fig = openfig([path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                                folder_output '/' predictions_folder '/' state_spaces{state_i} '/figures/Msfb_lengthall.fig']);
                            saveas(fig, [path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                                folder_output '/' predictions_folder '/' state_spaces{state_i} '/figures/Msfb_lengthall.jpeg']);
                            
                            %% length 3
                
                            % average across subjects/sessions
                            null_Msz_thr_pos = median(null_Msz_thr_pos_l3_allsub);
                            null_Msz_thr_neg = median(null_Msz_thr_neg_l3_allsub);
                            null_Msf_thr_pos = median(null_Msf_thr_pos_l3_allsub);
                            null_Msf_thr_neg = median(null_Msf_thr_neg_l3_allsub);
                            null_Msb_thr_pos = median(null_Msb_thr_pos_l3_allsub);
                            null_Msb_thr_neg = median(null_Msb_thr_neg_l3_allsub);

                            % avoid full overlap of threshold lines
                            if round(null_Msf_thr_pos,3) == round(null_Msb_thr_pos,3)
                                null_Msf_thr_pos = null_Msf_thr_pos + 0.005;
                            end

                            if round(null_Msf_thr_neg,3) == round(null_Msb_thr_neg,3)
                                null_Msf_thr_neg = null_Msf_thr_neg + 0.005;
                            end
            
            
                            Msz = median(Msz_l3_allsub, 2);
                            Msf = median(Msf_l3_allsub, 2);
                            Msb = median(Msb_l3_allsub, 2);

                            % plot Msf and Msb
                            figure;
                            plot((1:nlags)*5, Msf, 'Color', [0 0 0.8], 'LineWidth', 4); hold on;
                            plot((1:nlags)*5, Msb, 'Color', [0.8 0 0], 'LineWidth', 4); hold on;

                            yline(null_Msf_thr_pos,'--', 'Color', [0 0 0.8], 'LineWidth', 3);
                            yline(null_Msf_thr_neg,'--', 'Color', [0 0 0.8], 'LineWidth', 3);
                            yline(null_Msb_thr_pos,'--', 'Color', [0.8 0 0], 'LineWidth', 3);
                            yline(null_Msb_thr_neg,'--', 'Color', [0.8 0 0], 'LineWidth', 3);
%                             title(['Forward and Backward all lengths ' inputs{input_i} ' ' classifier ' ' state_spaces{state_i}])
%                             legend('forward', 'backward', 'Location', 'northwest');
%                             ylabel('sequenceness');
%                             xlabel('time lag (ms)');
                            yticks([-0.1 0 0.1 0.2]);
                            xticks([]);

                            % Increase the font size of the tick labels
                            ax = gca; % Get current axes                            
                            ax.FontSize = 20;
                
                            saveas(gcf, [path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                                folder_output '/' predictions_folder '/' state_spaces{state_i} '/figures/Msfb_length3.fig']);
                            fig = openfig([path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                                folder_output '/' predictions_folder '/' state_spaces{state_i} '/figures/Msfb_length3.fig']);
                            saveas(fig, [path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                                folder_output '/' predictions_folder '/' state_spaces{state_i} '/figures/Msfb_length3.jpeg']);
                            
                            %% length 4
                
                            % average across subjects/sessions
                            null_Msz_thr_pos = median(null_Msz_thr_pos_l4_allsub);
                            null_Msz_thr_neg = median(null_Msz_thr_neg_l4_allsub);
                            null_Msf_thr_pos = median(null_Msf_thr_pos_l4_allsub);
                            null_Msf_thr_neg = median(null_Msf_thr_neg_l4_allsub);
                            null_Msb_thr_pos = median(null_Msb_thr_pos_l4_allsub);
                            null_Msb_thr_neg = median(null_Msb_thr_neg_l4_allsub);
            
                            % avoid full overlap of threshold lines
                            if round(null_Msf_thr_pos,3) == round(null_Msb_thr_pos,3)
                                null_Msf_thr_pos = null_Msf_thr_pos + 0.005;
                            end

                            if round(null_Msf_thr_neg,3) == round(null_Msb_thr_neg,3)
                                null_Msf_thr_neg = null_Msf_thr_neg + 0.005;
                            end
            
                            
                            Msz = median(Msz_l4_allsub, 2);
                            Msf = median(Msf_l4_allsub, 2);
                            Msb = median(Msb_l4_allsub, 2);

                            % plot Msf and Msb
                            figure;
                            plot((1:nlags)*5, Msf, 'Color', [0 0 0.8], 'LineWidth', 4); hold on;
                            plot((1:nlags)*5, Msb, 'Color', [0.8 0 0], 'LineWidth', 4); hold on;

                            yline(null_Msf_thr_pos,'--', 'Color', [0 0 0.8], 'LineWidth', 3);
                            yline(null_Msf_thr_neg,'--', 'Color', [0 0 0.8], 'LineWidth', 3);
                            yline(null_Msb_thr_pos,'--', 'Color', [0.8 0 0], 'LineWidth', 3);
                            yline(null_Msb_thr_neg,'--', 'Color', [0.8 0 0], 'LineWidth', 3);
%                             title(['Forward and Backward all lengths ' inputs{input_i} ' ' classifier ' ' state_spaces{state_i}])
%                             legend('forward', 'backward', 'Location', 'northwest');
%                             ylabel('sequenceness');
%                             xlabel('time lag (ms)');
                            yticks([-0.1 0 0.1 0.2]);
                            xticks([]);

                            % Increase the font size of the tick labels
                            ax = gca; % Get current axes                            
                            ax.FontSize = 20;
                
                            saveas(gcf, [path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                                folder_output '/' predictions_folder '/' state_spaces{state_i} '/figures/Msfb_length4.fig']);
                            fig = openfig([path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                                folder_output '/' predictions_folder '/' state_spaces{state_i} '/figures/Msfb_length4.fig']);
                            saveas(fig, [path_results '/group_results/' inputs{input_i} '/' classifier '/' ...
                                folder_output '/' predictions_folder '/' state_spaces{state_i} '/figures/Msfb_length4.jpeg']);
                                            
                            close all

                        end
        
                    end
        
                end
        
            end
        
        end
    
    catch
    end
    
end

disp('Analysis done');

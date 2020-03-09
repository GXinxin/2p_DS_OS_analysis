clear; clc;

% path = 'E:\Lab\Data\2p\summaryData\newDataSet_042219\deep\';
% path = 'E:\Lab\Data\2p\summaryData\newDataSet_042219\deep\sameDeeper\';
% path = 'E:\Lab\Data\2p\tra2b\summary\';
% path = 'E:\Lab\Data\2p\FRMD7\newSummary_111419\';
% path = 'E:\Lab\Data\2p\FRMD7\P22\newSummary_111419\';
path = 'E:\Lab\Data\2p\summaryData\newSummary_111419\';
% path = 'E:\Lab\Data\2p\summaryData\p21\newSummaryData_111419\';

cd(path)

% tag = {'saline', 'gabazine_20x', 'gabazine_100x'};
tag = {'saline', 'gabazine', 'glutamate', 'nonInj'};
tag = tag([1, 4]);
% tag = {'ctrl', 'exp'};


thresh = [0.3 0.5]; % DSI th = 0.3 or th = 0.5
include_th = [0.2 0.3 0.4];
is_smooth = 0;

for depth_tag = 2 % 1 = superficial layer, 2 = deep layer
    if depth_tag == 1
        match_pattern = {'DS1', 'DS_1', 'DS1_1', 'DS_1_1', 'DS2', 'DS_2', 'DS2_1', 'DS_2_1'};
    else
        match_pattern = {'DS1_2', 'DS_1_2', 'DS2_2', 'DS_2_2'};
    end
%     param_tag = {'\p15_salineGabazine_', '*slc*'};
%     param_tag = {'\p15_salineGlutamate_', '*slc*'};
    param_tag = {'\p15_salineNoninj_', '*slc*'};
%     param_tag = {'\p21_frmd7_', '*rmd7*'};
%     param_tag = {'\p15_frmd7_', '*rmd7*'};
%     param_tag = {'\p21_salineGabazine_', '*slc*'};
%     param_tag = {'\p21_salineGlutamate_', '*slc*'};
%     param_tag = {'\p21_salineNoninj_', '*slc*'};
    colorVector = copper(4);
    colorVector = colorVector([1 4], :);
    
    
    for ii = 2 : length(include_th)
        clear is_include
        if is_smooth
            fd_name = [num2str(ii), '_smooth_depth', num2str(depth_tag)];
            mkdir(fd_name)
        else
            fd_name = [num2str(ii), '_noSmooth_depth', num2str(depth_tag)];
            mkdir(fd_name)
        end
        save_tag = [fd_name, param_tag{1}];
        
        % initialization
        DSI_total = cell(3, 2);
        OSI_total = cell(3, 2);
        DSI_fit_total = cell(3, 2);
        OSI_fit_total = cell(3, 2);
        DSI_global_total = cell(3, 2);
        OSI_global_total = cell(3, 2);
        time2peak_total = cell(3, 2);
        time2slope_total = cell(3, 2);
        cell_no = cell(3, 2);
        resp_Percent_all = cell(3, 2);
        percent_DS = cell(3, 2);
        percent_OS = cell(3, 2);
        percent_DS_fit = cell(3, 2);
        percent_OS_fit = cell(3, 2);
        DS_direction_total = cell(3, 2, 2);
        DS_direction_fit_total = cell(3, 2, 2);
        OS_direction_total = cell(3, 2, 2);
        OS_direction_fit_total = cell(3, 2, 2);
        
        
        
        
        for g = 1 : length(tag)
            
            cd([path, tag{g}])
            fdlist = dir(fullfile(param_tag{2}));
            
            count = 1;
            
            % for each animal
            for ff = 1 : length(fdlist)
                cd(fdlist(ff).name)
                sub_fdlist = dir(fullfile('*DS*'));
                
                subfd_id = [];
                for m = 1 : length(sub_fdlist)
                    for mm = 1 : length(match_pattern)
                        match_id = strfind(sub_fdlist(m).name, match_pattern{mm});
                        if ~isempty(match_id)
                            subfd_id = [subfd_id, m];
                        end
                    end
                end
                subfd_id = unique(subfd_id);
                
                
                % for each FOV per animal
                for f = 1 : length(subfd_id)
                    indexingM(:, count) = [ff; subfd_id(f)];
                    cd(sub_fdlist(subfd_id(f)).name)
                    if is_smooth
                        load('filtered_DSI_summary_window5.mat')
                    else
                        load('DSI_summary_window5.mat')
                    end
                    
                    % exclude FOVs with less than 20% responsive cells during
                    % the drift grating period
                    is_include{g}(ff, f) = (sum(is_respond(:, 2)) / length(is_respond(:, 2))) > include_th(ii);
                    if is_include{g}(ff, f)
                        for t = 1 : 3
                            cell_no{t, g}(:, count) = [sum(is_respond(:, t)); length(is_respond(:, t))];
                            resp_Percent_all{t, g}(count) = sum(is_respond(:, t)) / length(is_respond(:, t)); % percentage of cells that respond to each stimuli
                            percent_DS{t, g}(:, count) = DS_percent(:, t);
                            percent_OS{t, g}(:, count) = OS_percent(:, t);
                            percent_DS_fit{t, g}(:, count) = fit_DS_percent(:, t);
                            percent_OS_fit{t, g}(:, count) = fit_OS_percent(:, t);
                            
                            DSI_total{t, g} = [DSI_total{t, g}; DSI(id1{t}, t)];
                            OSI_total{t, g} = [OSI_total{t, g}; OSI(id1{t}, t)];
                            DSI_fit_total{t, g} = [DSI_fit_total{t, g}; DSI_fit(id2{t}, t)];
                            OSI_fit_total{t, g} = [OSI_fit_total{t, g}; OSI_fit(id2{t}, t)];
                            DSI_global_total{t, g} = [DSI_global_total{t, g}; DSI_global(id1{t}, t)];
                            OSI_global_total{t, g} = [OSI_global_total{t, g}; OSI_global(id1{t}, t)];
                            time2peak_total{t, g} = [time2peak_total{t, g}, time2peak{t}];
                            time2slope_total{t, g} = [time2slope_total{t, g}, time2slope{t}];
                            
                            DSI_individual{t, g, count} = DSI(id1{t}, t);
                            OSI_individual{t, g, count} = OSI(id1{t}, t);
                            DSI_fit_individual{t, g, count} = DSI_fit(id2{t}, t);
                            OSI_fit_individual{t, g, count} = OSI_fit(id2{t}, t);
                            DSI_global_individual{t, g, count} = DSI_global(id1{t}, t);
                            OSI_global_individual{t, g, count} = OSI_global(id1{t}, t);
                            time2peak_individual{t, g, count} = time2peak{t};
                            time2slope_individual{t, g, count} = time2slope{t};
                            
                            for th = 1 : 2
                                DS_direction_total{t, th, g} = [DS_direction_total{t, th, g}; DS_direction{t}{th}];
                                DS_direction_fit_total{t, th, g} = [DS_direction_fit_total{t, th, g}, fit_DS_direction{t}{th}];
                                OS_direction_total{t, th, g} = [OS_direction_total{t, th, g}; OS_direction{t}{th}];
                                OS_direction_fit_total{t, th, g} = [OS_direction_fit_total{t, th, g}, fit_OS_direction{t}{th}];
                                
                                % also saves DS and OS properties for individual FOVs
                                DS_direction_individual{t, th, g}{count} = DS_direction{t}{th};
                                DS_direction_individual_fit{t, th, g}{count} = fit_DS_direction{t}{th};
                                OS_direction_individual{t, th, g}{count} = OS_direction{t}{th};
                                OS_direction_individual_fit{t, th, g}{count} = fit_OS_direction{t}{th};
                            end
                        end
                        median_individual{g}(:, :, count) = median_values; % col1: DSI, col2: DSI_fit, col3: DSI_global, col4: OSI, col5: OSI_fit, col6: OSI_global
                        count = count + 1;
                    end
                    cd ..
                end
                cd ..
            end
            
            
            for t = 1 : 3
                resp_percent_total{t}(g) = sum(cell_no{t, g}(1, :)) / sum(cell_no{t, g}(2, :));
                
                median_total{g}(t, 1) = median(DSI_total{t, g});
                median_total{g}(t, 2) = median(DSI_fit_total{t, g});
                median_total{g}(t, 3) = median(DSI_global_total{t, g});
                median_total{g}(t, 4) = median(OSI_total{t, g});
                median_total{g}(t, 5) = median(OSI_fit_total{t, g});
                median_total{g}(t, 6) = median(OSI_global_total{t, g});
                median_total{g}(t, 7) = median(time2peak_total{t, g});
                median_total{g}(t, 8) = median(time2slope_total{t, g});
                
                for th = 1 : 2
                    DS_percent_total{g}(th, t) = sum(DSI_total{t, g} > thresh(th));
                    DS_percent_fit_total{g}(th, t) = sum(DSI_fit_total{t, g} > thresh(th));
                    OS_percent_total{g}(th, t) = sum(OSI_total{t, g} > thresh(th));
                    OS_percent_fit_total{g}(th, t) = sum(OSI_fit_total{t, g} > thresh(th));
                end
                
            end
            
            cd ..
            
        end
        
        
        
        
        %% DSI, OSI, time2peak/time2slope
        % titleText = {'glutamate blocker', 'cholinergic blocker', 'saline'};
        % titleText = {'glutamate blocker', 'saline'};
        titleText = {'ctrl', 'exp'};
        
        
        for cc = 1 : 8
            
            h = figure; set(h, 'position', [0 0 1000, 1000])
            count = 1;
            for t = 1 : 3
                for g = 1:length(tag)
                    
                    switch cc
                        case 1
                            data = DSI_total{t, g};
                            cc_tag = 'DSI';
                        case 2
                            data = DSI_fit_total{t, g};
                            median_data = median_total{g}(t, 1);
                            cc_tag = 'DSI_fit';
                        case 3
                            data = DSI_global_total{t, g};
                            cc_tag = 'DSI_global';
                        case 4
                            data = OSI_total{t, g};
                            cc_tag = 'OSI';
                        case 5
                            data = OSI_fit_total{t, g};
                            cc_tag = 'OSI_fit';
                        case 6
                            data = OSI_global_total{t, g};
                            cc_tag = 'OSI_global';
                        case 7
                            data = time2peak_total{t, g};
                            cc_tag = 'time2peak';
                        case 8
                            data = time2slope_total{t, g};
                            cc_tag = 'time2slope';
                    end
                    median_data = median_total{g}(t, cc);
                    if g == 1
                        pre_data = data;
                    end
                    
                    
                    subplot(3, 2, count)
                    if cc > 6
                        [counts(g, :), centers(g, :)] = hist(data, 20);
                    else
                        [counts(g, :), centers(g, :)] = hist(data, 0.025:0.05:1);
                    end
                    proportion(g, :) = counts(g, :) / length(data);
                    bar(centers(g, :), proportion(g, :), 'barWidth', 0.9, 'FaceColor', colorVector(g, :))
                    line(median_data * [1 1], [0 0.3], 'color', 'r', 'lineWidth', 2)
                    if cc > 6
                        xlim([0 150])
                    else
                        xlim([0 1])
                    end
                    ylim([0 0.45])
                    ylabel('proportion', 'FontSize', 12, 'FontWeight', 'Bold')
                    xlabel(cc_tag, 'FontSize', 12, 'FontWeight', 'Bold')
                    title([titleText{g}, ' median=', num2str(median_data)], 'FontSize', 15)
                    count = count + 1;
                    box off
                end
                [p1(t, cc), hh1(t, cc)] = ttest2(pre_data, data);
            end
            saveas(h, [save_tag, [cc_tag, '_hist.png']])
            
        end
        
        
        
        
        %% plots per FOV or per animal
        plot_FOV_summary_oldVersion;
        plot_animal_summary_oldVersion;
        
        
    end
end

%%
% % preferred direction distribution for each animals, separated by
% % conditions
% for g = 1 : 2
%     h = figure;
%     set(h, 'position', [0 0 900 900])
%     for ff = 1 : length(DS_direction_individual{g})
%         [counts_d_ind(ff, :), centers_d_ind(ff, :)] = hist(DS_direction_individual{g}{ff}, 1:8);
%         proportion_d_ind(ff, :) = counts_d_ind(ff, :) / length(DS_direction_individual{g}{ff});
%         polarplot(((0 : 1/4 : 2) - 1/2) * pi, [proportion_d_ind(ff, :), proportion_d_ind(ff, 1)], 'Color', colorVector(5-g, :), ...
%             'linewidth', 1)
%         hold on
%     end
%     [counts_d(g, :), centers_d(g, :)] = hist(DS_direction_total{g}, 1:8);
%     proportion_d(g, :) = counts_d(g, :) / length(DS_direction_total{g});
%     polarplot(((0 : 1/4 : 2) - 1/2) * pi, [proportion_d(g, :), proportion_d(g, 1)], 'Color', 'k', ...
%         'linewidth', 4)
%     title(tag{g})
%     rlim([0, 0.7])
%     set(gca,'linew',2, 'FontSize', 22)
%     saveas(h, [save_tag, 'distribution_', tag{g}, '.png'])
%
%     circ_rtest(DS_direction_total{g}*pi/4)
% end
%








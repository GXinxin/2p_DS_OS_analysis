% combine FOVs from one animal
for g = 1:length(tag)
    a = 1;
    fov_id = 1;
    for ff = 1 : length(is_include{g})
        if sum(is_include{g}(ff, :)) > 0
            
            if sum(is_include{g}(ff, :)) == 2
                select_fov = [fov_id, fov_id+1]
                fov_id = fov_id + 2;
            else
                select_fov = fov_id
                fov_id = fov_id + 1;
            end
            
            for t = 1 : 3
                resp_Percent_animal{t, g}(a) = sum(cell_no{t, g}(1, select_fov))/sum(cell_no{t, g}(2, select_fov));
                percent_DS_animal{t, g}(:, a) = sum(percent_DS{t, g}(:, select_fov) .* repmat(cell_no{t, g}(1, select_fov), 2, 1), 2) ...
                    ./ repmat(sum(cell_no{t, g}(1, select_fov)), 2, 1);
                percent_DS_fit_animal{t, g}(:, a) = sum(percent_DS_fit{t, g}(:, select_fov) .* repmat(cell_no{t, g}(1, select_fov), 2, 1), 2) ...
                    ./ repmat(sum(cell_no{t, g}(1, select_fov)), 2, 1);
                percent_OS_animal{t, g}(:, a) = sum(percent_OS{t, g}(:, select_fov) .* repmat(cell_no{t, g}(1, select_fov), 2, 1), 2) ...
                    ./ repmat(sum(cell_no{t, g}(1, select_fov)), 2, 1);
                percent_OS_fit_animal{t, g}(:, a) = sum(percent_OS_fit{t, g}(:, select_fov) .* repmat(cell_no{t, g}(1, select_fov), 2, 1), 2) ...
                    ./ repmat(sum(cell_no{t, g}(1, select_fov)), 2, 1);
%                 percent_OS_DS_animal{t, g}(:, a) = sum(percent_OS_DS{t, g}(:, select_fov) .* repmat(cell_no{t, g}(1, select_fov), 2, 1), 2) ...
%                     ./ repmat(sum(cell_no{t, g}(1, select_fov)), 2, 1);
%                 percent_OS_only_animal{t, g}(:, a) = sum(percent_OS_only{t, g}(:, select_fov) .* repmat(cell_no{t, g}(1, select_fov), 2, 1), 2) ...
%                     ./ repmat(sum(cell_no{t, g}(1, select_fov)), 2, 1);
%                 percent_OS_DS_fit_animal{t, g}(:, a) = sum(percent_OS_DS_fit{t, g}(:, select_fov) .* repmat(cell_no{t, g}(1, select_fov), 2, 1), 2) ...
%                     ./ repmat(sum(cell_no{t, g}(1, select_fov)), 2, 1);
%                 percent_OS_only_fit_animal{t, g}(:, a) = sum(percent_OS_only_fit{t, g}(:, select_fov) .* repmat(cell_no{t, g}(1, select_fov), 2, 1), 2) ...
%                     ./ repmat(sum(cell_no{t, g}(1, select_fov)), 2, 1);
%                 
                
                
                
                
                if length(select_fov) > 1
                    median_individual_animal{g}(t, 1, a) = median([DSI_individual{t, g, select_fov(1)}; DSI_individual{t, g, select_fov(2)}]);
                    median_individual_animal{g}(t, 2, a) = median([DSI_fit_individual{t, g, select_fov(1)}; DSI_fit_individual{t, g, select_fov(2)}]);
                    median_individual_animal{g}(t, 3, a) = median([DSI_global_individual{t, g, select_fov(1)}; DSI_global_individual{t, g, select_fov(2)}]);
                    median_individual_animal{g}(t, 4, a) = median([OSI_individual{t, g, select_fov(1)}; OSI_individual{t, g, select_fov(2)}]);
                    median_individual_animal{g}(t, 5, a) = median([OSI_fit_individual{t, g, select_fov(1)}; OSI_fit_individual{t, g, select_fov(2)}]);
                    median_individual_animal{g}(t, 6, a) = median([OSI_global_individual{t, g, select_fov(1)}; OSI_global_individual{t, g, select_fov(2)}]);
                    median_individual_animal{g}(t, 7, a) = median([time2peak_individual{t, g, select_fov(1)}, time2peak_individual{t, g, select_fov(2)}]);
                    median_individual_animal{g}(t, 8, a) = median([time2slope_individual{t, g, select_fov(1)}, time2slope_individual{t, g, select_fov(2)}]);
                    
                    
                    if t == 2
%                         for th = 1 : 2
%                             pair_counts_OS_animal{th, g}{a} = {};
%                             pair_counts_fit_OS_animal{th, g}{a} = {};
%                             pair_counts_DS_animal{th, g}{a} = {};
%                             pair_counts_fit_DS_animal{th, g}{a} = {};
%                             pair_counts_OS_animal_new{th, g}{a} = {};
%                             pair_counts_fit_OS_animal_new{th, g}{a} = {};
%                             pair_counts_DS_animal_new{th, g}{a} = {};
%                             pair_counts_fit_DS_animal_new{th, g}{a} = {};
%                             
%                             
%                             for range = 1 : 6
%                                 tmp = pair_counts_OS{t, th, g};
%                                 pair_counts_OS_animal{th, g}{a}{range} = combineFOV(tmp{select_fov(1)}{range}, tmp{select_fov(2)}{range});
%                                 pair_counts_fit_OS_animal{th, g}{a}{range} = combineFOV(pair_counts_fit_OS{t, th, g}{select_fov(1)}{range}, pair_counts_fit_OS{t, th, g}{select_fov(2)}{range});
%                                 pair_counts_DS_animal{th, g}{a}{range} = combineFOV(pair_counts_DS{t, th, g}{select_fov(1)}{range}, pair_counts_DS{t, th, g}{select_fov(2)}{range});
%                                 pair_counts_fit_DS_animal{th, g}{a}{range} = combineFOV(pair_counts_fit_DS{t, th, g}{select_fov(1)}{range}, pair_counts_fit_DS{t, th, g}{select_fov(2)}{range});
%                             end
%                             
%                             for range = 1 : 12
%                                 pair_counts_OS_animal_new{th, g}{a}{range} = combineFOV(pair_counts_OS_new{t, th, g}{select_fov(1)}{range}, pair_counts_OS_new{t, th, g}{select_fov(2)}{range});
%                                 pair_counts_fit_OS_animal_new{th, g}{a}{range} = combineFOV(pair_counts_fit_OS_new{t, th, g}{select_fov(1)}{range}, pair_counts_fit_OS_new{t, th, g}{select_fov(2)}{range});
%                                 pair_counts_DS_animal_new{th, g}{a}{range} = combineFOV(pair_counts_DS_new{t, th, g}{select_fov(1)}{range}, pair_counts_DS_new{t, th, g}{select_fov(2)}{range});
%                                 pair_counts_fit_DS_animal_new{th, g}{a}{range} = combineFOV(pair_counts_fit_DS_new{t, th, g}{select_fov(1)}{range}, pair_counts_fit_DS_new{t, th, g}{select_fov(2)}{range});
%                             end
%                         end
                    end
                else
                    median_individual_animal{g}(t, 1, a) = median(DSI_individual{t, g, select_fov});
                    median_individual_animal{g}(t, 2, a) = median(DSI_fit_individual{t, g, select_fov});
                    median_individual_animal{g}(t, 3, a) = median(DSI_global_individual{t, g, select_fov});
                    median_individual_animal{g}(t, 4, a) = median(OSI_individual{t, g, select_fov});
                    median_individual_animal{g}(t, 5, a) = median(OSI_fit_individual{t, g, select_fov});
                    median_individual_animal{g}(t, 6, a) = median(OSI_global_individual{t, g, select_fov});
                    median_individual_animal{g}(t, 7, a) = median(time2peak_individual{t, g, select_fov});
                    median_individual_animal{g}(t, 8, a) = median(time2slope_individual{t, g, select_fov});
                    
                    if t == 2
                        for th = 1 : 2
%                             pair_counts_OS_animal{th, g}{a} = pair_counts_OS{t, th, g}{select_fov} ;
%                             pair_counts_fit_OS_animal{th, g}{a} = pair_counts_fit_OS{t, th, g}{select_fov};
%                             pair_counts_DS_animal{th, g}{a} = pair_counts_DS{t, th, g}{select_fov};
%                             pair_counts_fit_DS_animal{th, g}{a} = pair_counts_fit_DS{t, th, g}{select_fov};
%                             
                            
%                             pair_counts_OS_animal_new{th, g}{a} = pair_counts_OS_new{t, th, g}{select_fov} ;
%                             pair_counts_fit_OS_animal_new{th, g}{a} = pair_counts_fit_OS_new{t, th, g}{select_fov};
%                             pair_counts_DS_animal_new{th, g}{a} = pair_counts_DS_new{t, th, g}{select_fov};
%                             pair_counts_fit_DS_animal_new{th, g}{a} = pair_counts_fit_DS_new{t, th, g}{select_fov};
                        end
                    end
                    
                end
                
                
            end
            a = a + 1;
        end
    end
end


%% responding cell percentage
isSingleTh = 0;
range = [1, 1, 1, 1, 1, 1, 1, 1, 1];

for cc = 1 : 5
    h = figure;
    set(h, 'position', [0 0 1500 450])
    
    for t = 1 : 3
        subplot(1, 3, t)
        
        for g = 1:length(tag)
            switch cc
                case 1
                    data = resp_Percent_animal{t, g};
                    cc_tag = 'responding';
                case 2
                    data = percent_DS_animal{t, g};
                    cc_tag = 'DSI';
                case 3
                    data = percent_DS_fit_animal{t, g};
                    cc_tag = 'DSI_fit';
                case 4
                    data = percent_OS_animal{t, g};
                    cc_tag = 'OSI';
                case 5
                    data = percent_OS_fit_animal{t, g};
                    cc_tag = 'OSI_fit';
                case 6
                    data = percent_OS_DS_animal{t, g};
                    cc_tag = 'OS_DS';
                case 7
                    data = percent_OS_only_animal{t, g};
                    cc_tag = 'OS_only';
                case 8
                    data = percent_OS_DS_fit_animal{t, g};
                    cc_tag = 'OS_DS_fit';
                case 9
                    data = percent_OS_only_fit_animal{t, g};
                    cc_tag = 'OS_only_fit';
            end
            
            if g == 1
                pre_data = data;
            end
            
            if cc == 1
                scatter(g * ones(1, length(data)), data, 'MarkerEdgeColor', colorVector(g, :)', 'LineWidth', 2)
                hold on
            else
                if isSingleTh
                    scatter(g * ones(1, length(data(2, :))), data(2, :), 'MarkerEdgeColor', colorVector(g, :)', 'LineWidth', 2)
                    hold on
                else
                    scatter(g * ones(1, length(data(1, :))), data(1, :), 'MarkerEdgeColor', colorVector(g, :)', 'LineWidth', 2)
                    hold on
                    scatter((g+2) * ones(1, length(data(2, :))), data(2, :), 'MarkerEdgeColor', colorVector(g, :)', 'LineWidth', 2)
                    hold on
                end
            end
        end
        
        clear tmp
        if cc == 1
            tmp{1} = pre_data;
            tmp{2} = data;
            colors = colorVector([1, 2], :);
            ftag = tag;
        else
            tmp{1} = pre_data(1, :);
            tmp{3} = pre_data(2, :);
            tmp{2} = data(1, :);
            tmp{4} = data(2, :);
            colors = colorVector([1, 2, 1, 2], :);
            ftag = {[tag{1}, '1'], [tag{2}, '1'], [tag{1}, '2'], [tag{2}, '2']};
        end
        
        col=@(x)reshape(x,numel(x),1);
        boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
        if isSingleTh
            if cc > 1
                boxplot2(tmp(3:4), 'Width', .6, 'Colors', colors(1:2, :), 'Labels', tag, 'whisker', 1000)
            else
                boxplot2(tmp, 'Width', .6, 'Colors', colors, 'Labels', ftag, 'whisker', 1000)
            end
        else
            boxplot2(tmp, 'Width', .6, 'Colors', colors, 'Labels', ftag, 'whisker', 1000)
        end
        
        title([cc_tag, ' t', num2str(t)])
        ylabel('proportion')
        set(findobj(gca,'type','line'),'linew',3)
        set(gca,'linew',4, 'FontSize', 14)
        if cc == 1
            xlim([0 length(tag)+1])
        else
            if isSingleTh
                xlim([0 length(tag)+1])
            else
                xlim([0 2*length(tag)+1])
            end
        end
        ylim([0 range(cc)])
        box off
        saveas(h, [save_tag, cc_tag, '_percent_perAnimal_onlyTh2.png'])
        
        
        if cc == 1
            [p2{1}(t, cc), hh2{1}(t, cc)] = ttest2(tmp{1}, tmp{2});
            p2{2}(t, cc) = NaN;
            hh2{2}(t, cc) = NaN;
        else
            [p2{1}(t, cc), hh2{1}(t, cc)] = ttest2(tmp{1}, tmp{2});
            [p2{2}(t, cc), hh2{2}(t, cc)] = ttest2(tmp{3}, tmp{4});
        end
    end
end




%% individual median values of DSI/OSI
% col1: DSI, col2: DSI_fit, col3: DSI_global, col4: OSI, col5: OSI_fit, col6: OSI_global
range = [0.6, 0.6, 0.35, 1, 1, 0.5];
for cc = 1 : 6
    
    colors = colorVector([1, 2], :);
    ftag = tag;
    
    h = figure;
    set(h, 'position', [0 0 1200 450])
    for t = 1 : 3
        subplot(1, 3, t)
        for g = 1:length(tag)
            switch cc
                case 1
                    cc_tag = 'DSI';
                case 2
                    cc_tag = 'DSI_fit';
                case 3
                    cc_tag = 'DSI_global';
                case 4
                    cc_tag = 'OSI';
                case 5
                    cc_tag = 'OSI_fit';
                case 6
                    cc_tag = 'OSI_global';
            end
            
            data = squeeze(median_individual_animal{g}(t, cc, :));
            scatter(g * ones(1, length(data)), data, 'MarkerEdgeColor', colorVector(g, :)', 'LineWidth', 2)
            hold on
        end
        
        clear tmp
        tmp{1} = squeeze(median_individual_animal{1}(t, cc, :));
        tmp{2} = squeeze(median_individual_animal{2}(t, cc, :));
        col=@(x)reshape(x,numel(x),1);
        boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
        boxplot2(tmp, 'Width', .6, 'Colors', colors, 'Labels', ftag, 'whisker', 1000)
        box off
        
        title([cc_tag, ' t', num2str(t)])
        xlim([0 length(tag)+1]); ylim([0 range(cc)])
        ylabel('proportion')
        set(findobj(gca,'type','line'),'linew',3)
        set(gca,'linew',4, 'FontSize', 14)
        
        [p3(t, cc), hh3(t, cc)] = ttest2(tmp{1}, tmp{2});
        
    end
    saveas(h, [save_tag, cc_tag, '_median_animal.png'])
    
end



%% spatial regularity
% colors = colorVector([1, 2], :);
% title_tag = {'OS', 'OS(f)', 'DS', 'DS(f)'};
% for t = 2
%     for cc = 1 : 4
%         
%         clear pair_mean pair_sem
%         
%         switch cc
%             case 1
%                 pair_c = pair_counts_OS_animal;
%             case 2
%                 pair_c = pair_counts_fit_OS_animal;
%             case 3
%                 pair_c = pair_counts_DS_animal;
%             case 4
%                 pair_c = pair_counts_fit_DS_animal;
%         end
%         
%         
%         for th = 1 : 2
%             total_pair = length(dist_p{cc, t, th});
%             
%             % plot difference vs proportion of pairs in 6 separate distance bins
%             h = figure; set(h, 'position', [0 0 1000 1000])
%             for dd = 1 : 6
%                 pair_mat{1} = []; pair_mat{2} = [];
%                 for fov_n = 1 : length(pair_c{th, 1})
%                     tmp = pair_c{th, 1};
%                     if ~isempty(tmp{fov_n}) && sum(tmp{fov_n}{dd}) > 0
%                         pair_mat{1} = [pair_mat{1}; tmp{fov_n}{dd}/sum(tmp{fov_n}{dd})];
%                     end
%                 end
%                 
%                 for fov_n = 1 : length(pair_c{th, 2})
%                     tmp = pair_c{th, 2};
%                     if ~isempty(tmp{fov_n}) && sum(tmp{fov_n}{dd}) > 0
%                         pair_mat{2} = [pair_mat{2}; tmp{fov_n}{dd}/sum(tmp{fov_n}{dd})];
%                     end
%                 end
%                 pair_mean(1, :) = mean(pair_mat{1});
%                 pair_sem(1, :) = std(pair_mat{1})/sqrt(size(pair_mat{1}, 2));
%                 pair_mean(2, :) = mean(pair_mat{2});
%                 pair_sem(2, :) = std(pair_mat{2})/sqrt(size(pair_mat{2}, 2));
%                 
%                 for i = 1 : size(pair_mat{1}, 2)
%                     if ~isempty(pair_mat{2}) && ~isempty(pair_mat{1})
%                         [p4{cc, t, th}(dd, i), hh4{cc, t, th}(dd, i)] = ttest2(pair_mat{1}(:, i), pair_mat{2}(:, i));
%                     end
%                 end
%                 
%                 
%                 subplot(2, 3, dd)
%                 errorbar(1:size(pair_mean, 2), pair_mean(1, :), pair_sem(1, :), 'color', colors(1, :), 'linewidth', 1); hold on
%                 errorbar(1:size(pair_mean, 2), pair_mean(2, :), pair_sem(2, :), 'color', colors(2, :), 'linewidth', 1);
%                 title(['d', num2str(dd), ' ', title_tag{cc}])
%                 box off
%             end
%             saveas(h, [save_tag, 'Spatial_regularity_t', num2str(t), '_', title_tag{cc}, num2str(th), '_animal.png'])
%         end
%     end
% end



%% plot distance vs difference
% t = 2;
% for cc = 1 : 4
%     
%     clear pair_mean pair_sem
%     
%     switch cc
%         case 1
%             pair_c = pair_counts_OS_animal_new;
%         case 2
%             pair_c = pair_counts_fit_OS_animal_new;
%         case 3
%             pair_c = pair_counts_DS_animal_new;
%         case 4
%             pair_c = pair_counts_fit_DS_animal_new;
%     end
%     
%     
%     for th = 1:2
%         h = figure;
%         
%         for g = 1 : 2
%             tmp = pair_c{th, g};
%             for a = 1 : length(tmp)
%                 individual_pair = tmp{a};
%                 for dd = 1 : 10
%                     if ~isempty(individual_pair{dd})
%                         dist_diff{g}(dd, a) = (individual_pair{dd} * (1 : length(individual_pair{dd}))') / sum(individual_pair{dd});
%                     end
%                 end
%             end
%             
%             mean_dist_pair(:, g) = mean(dist_diff{g}, 2);
%             std_dist_pair(:, g) = std(dist_diff{g}, [], 2);
%         end
%         
%         for dd = 1 : 10
%             [p5{th}(cc, dd), hh5{th}(cc, dd)] = ttest2(dist_diff{1}(dd, :), dist_diff{2}(dd, :));
%         end
%         
%         errorbar(1:10, mean_dist_pair(:, 1), std_dist_pair(:, 1)/sqrt(size(dist_diff{1}, 2)), 'color', colorVector(1, :), 'lineWidth', 2); hold on
%         errorbar(1:10, mean_dist_pair(:, 2), std_dist_pair(:, 2)/sqrt(size(dist_diff{2}, 2)), 'color', colorVector(2, :), 'lineWidth', 2); hold on
%         title(title_tag{cc})
%         saveas(h, [save_tag, 'Spatial_distVsSimilarity_t', num2str(t), '_', title_tag{cc}, num2str(th), '_animal.png'])
%     end
% end






%%
close all
fn_tmp = {'noSmoothing', 'smoothing'};
save([save_tag, 'stats_', fn_tmp{is_smooth+1}, '_th', num2str(ii), '_animal.mat'], 'hh1', 'hh2', 'hh3', 'is_include')




% this function combines spatial data from two FOVs of the same animal
function combined_data = combineFOV(data1, data2)
if ~isempty(data1) && ~isempty(data2)
    combined_data = data1 + data2;
else
    if ~isempty(data1)
        combined_data = data1;
    elseif ~isempty(data2)
        combined_data = data2;
    else
        combined_data = [];
    end
end
end

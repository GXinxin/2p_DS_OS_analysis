%% compute DSI, DSI_global, OSI and OSI_global for cells preferring different directions without thresholding
% threshold = [0.1:0.1:0.6, 0.15, 0.25];
threshold = 0:0.05:0.5;
t = 2;
for g = 1 : length(tag)
    % directions
    fov_n = size(DSI_fit_perDir{g}{1});
    for d = 1 : 8
        for ii = 1 : 3
            total_index_perDir{ii}{g}{d} = [];
        end
        for f = 1 : fov_n(2)
            total_index_perDir{1}{g}{d} = [total_index_perDir{1}{g}{d}; DSI_perDir{g}{d}{f}];
            total_index_perDir{2}{g}{d} = [total_index_perDir{2}{g}{d}; DSI_fit_perDir{g}{d}{f}];
            total_index_perDir{3}{g}{d} = [total_index_perDir{3}{g}{d}; DSI_global_perDir{g}{d}{f}];
            median_index_perDir{1}{g}(f, d) = median(DSI_perDir{g}{d}{f});
            median_index_perDir{2}{g}(f, d) = median(DSI_fit_perDir{g}{d}{f});
            median_index_perDir{3}{g}(f, d) = median(DSI_global_perDir{g}{d}{f});
            for tt = 1 : length(threshold)
                percent_perDir{1}{tt}{g}(f, d, 1) = sum(DSI_perDir{g}{d}{f} > threshold(tt))/length(DSI_perDir{g}{d}{f});
                percent_perDir{1}{tt}{g}(f, d, 2) = sum(DSI_perDir{g}{d}{f} > threshold(tt))/cell_no{g}(f);
                percent_perDir{2}{tt}{g}(f, d, 1) = sum(DSI_fit_perDir{g}{d}{f} > threshold(tt))/length(DSI_fit_perDir{g}{d}{f});
                percent_perDir{2}{tt}{g}(f, d, 2) = sum(DSI_fit_perDir{g}{d}{f} > threshold(tt))/cell_no{g}(f);
                percent_perDir{3}{tt}{g}(f, d, 1) = sum(DSI_global_perDir{g}{d}{f} > threshold(tt))/length(DSI_global_perDir{g}{d}{f});
                percent_perDir{3}{tt}{g}(f, d, 2) = sum(DSI_global_perDir{g}{d}{f} > threshold(tt))/cell_no{g}(f);
                DScell_no_thresh{1}{tt}{g}(f, d) = sum(DSI_perDir{g}{d}{f} > threshold(tt));
                DScell_no_thresh{2}{tt}{g}(f, d) = sum(DSI_fit_perDir{g}{d}{f} > threshold(tt));
                DScell_no_thresh{3}{tt}{g}(f, d) = sum(DSI_global_perDir{g}{d}{f} > threshold(tt));
            end
        end
        
        fov_id = 1;
        %         count = 1;
        for ff = 1 : length(is_include{g})
            if sum(is_include{g}(ff, :)) > 0
                
                if sum(is_include{g}(ff, :)) == 2
                    select_fov = [fov_id, fov_id+1];
                    fov_id = fov_id + 2;
                else
                    select_fov = fov_id;
                    fov_id = fov_id + 1;
                end
                median_index_perDir_animal{1}{g}(ff, d) = median(cell2mat(DSI_perDir{g}{d}(select_fov)'));
                median_index_perDir_animal{2}{g}(ff, d) = median(cell2mat(DSI_fit_perDir{g}{d}(select_fov)'));
                median_index_perDir_animal{3}{g}(ff, d) = median(cell2mat(DSI_global_perDir{g}{d}(select_fov)'));
                
                
                for tt = 1 : length(threshold)
                    if length(select_fov) < 2
                        tmp = cell2mat(DSI_perDir{g}{d}(select_fov)');
                        percent_perDir_animal{1}{tt}{g}(ff, d, 1) = sum(tmp > threshold(tt))/length(tmp);
                        percent_perDir_animal{1}{tt}{g}(ff, d, 2) = sum(tmp > threshold(tt))/(cell_no{t, g}(3, select_fov(1)));
                        tmp = cell2mat(DSI_fit_perDir{g}{d}(select_fov)');
                        percent_perDir_animal{2}{tt}{g}(ff, d, 1) = sum(tmp > threshold(tt))/length(tmp);
                        percent_perDir_animal{2}{tt}{g}(ff, d, 2) = sum(tmp > threshold(tt))/(cell_no{t, g}(3, select_fov(1)));
                        tmp = cell2mat(DSI_global_perDir{g}{d}(select_fov)');
                        percent_perDir_animal{3}{tt}{g}(ff, d, 1) = sum(tmp > threshold(tt))/length(tmp);
                        percent_perDir_animal{3}{tt}{g}(ff, d, 2) = sum(tmp > threshold(tt))/(cell_no{t, g}(3, select_fov(1)));
                    else
                        tmp = cell2mat(DSI_perDir{g}{d}(select_fov)');
                        percent_perDir_animal{1}{tt}{g}(ff, d, 1) = sum(tmp > threshold(tt))/length(tmp);
                        percent_perDir_animal{1}{tt}{g}(ff, d, 2) = sum(tmp > threshold(tt))/(cell_no{t, g}(3, select_fov(1)) + cell_no{t, g}(3, select_fov(2)));
                        tmp = cell2mat(DSI_fit_perDir{g}{d}(select_fov)');
                        percent_perDir_animal{2}{tt}{g}(ff, d, 1) = sum(tmp > threshold(tt))/length(tmp);
                        percent_perDir_animal{2}{tt}{g}(ff, d, 2) = sum(tmp > threshold(tt))/(cell_no{t, g}(3, select_fov(1)) + cell_no{t, g}(3, select_fov(2)));
                        tmp = cell2mat(DSI_global_perDir{g}{d}(select_fov)');
                        percent_perDir_animal{3}{tt}{g}(ff, d, 1) = sum(tmp > threshold(tt))/length(tmp);
                        percent_perDir_animal{3}{tt}{g}(ff, d, 2) = sum(tmp > threshold(tt))/(cell_no{t, g}(3, select_fov(1)) + cell_no{t, g}(3, select_fov(2)));
                    end
                end
                %                 count = count + 1;
            end
            
        end
    end
    
    
    % orientations
    for o = 1 : 4
        for ii = 1 : 3
            total_index_perOri{ii}{g}{o} = [];
        end
        for f = 1 : fov_n(2)
            total_index_perOri{1}{g}{o} = [total_index_perOri{1}{g}{o}; OSI_perOri{g}{o}{f}];
            total_index_perOri{2}{g}{o} = [total_index_perOri{2}{g}{o}; OSI_fit_perOri{g}{o}{f}];
            total_index_perOri{3}{g}{o} = [total_index_perOri{3}{g}{o}; OSI_global_perOri{g}{o}{f}];
            median_index_perOri{1}{g}(f, o) = median(OSI_perOri{g}{o}{f});
            median_index_perOri{2}{g}(f, o) = median(OSI_fit_perOri{g}{o}{f});
            median_index_perOri{3}{g}(f, o) = median(OSI_global_perOri{g}{o}{f});
            for tt = 1 : length(threshold)
                percent_perOri{1}{tt}{g}(f, o, 1) = sum(OSI_perOri{g}{o}{f} > threshold(tt))/length(OSI_perOri{g}{o}{f});
                percent_perOri{1}{tt}{g}(f, o, 2) = sum(OSI_perOri{g}{o}{f} > threshold(tt))/cell_no{g}(f);
                percent_perOri{2}{tt}{g}(f, o, 1) = sum(OSI_fit_perOri{g}{o}{f} > threshold(tt))/length(OSI_fit_perOri{g}{o}{f});
                percent_perOri{2}{tt}{g}(f, o, 2) = sum(OSI_fit_perOri{g}{o}{f} > threshold(tt))/cell_no{g}(f);
                percent_perOri{3}{tt}{g}(f, o, 1) = sum(OSI_global_perOri{g}{o}{f} > threshold(tt))/length(OSI_global_perOri{g}{o}{f});
                percent_perOri{3}{tt}{g}(f, o, 2) = sum(OSI_global_perOri{g}{o}{f} > threshold(tt))/cell_no{g}(f);
                OScell_no_thresh{1}{tt}{g}(f, o) = sum(OSI_perOri{g}{o}{f} > threshold(tt));
                OScell_no_thresh{2}{tt}{g}(f, o) = sum(OSI_fit_perOri{g}{o}{f} > threshold(tt));
                OScell_no_thresh{3}{tt}{g}(f, o) = sum(OSI_global_perOri{g}{o}{f} > threshold(tt));
            end
        end
        
        
        fov_id = 1;
        %         count = 1;
        for ff = 1 : length(is_include{g})
            if sum(is_include{g}(ff, :)) > 0
                
                if sum(is_include{g}(ff, :)) == 2
                    select_fov = [fov_id, fov_id+1];
                    fov_id = fov_id + 2;
                else
                    select_fov = fov_id;
                    fov_id = fov_id + 1;
                end
                median_index_perOri_animal{1}{g}(ff, o) = median(cell2mat(OSI_perOri{g}{o}(select_fov)'));
                median_index_perOri_animal{2}{g}(ff, o) = median(cell2mat(OSI_fit_perOri{g}{o}(select_fov)'));
                median_index_perOri_animal{3}{g}(ff, o) = median(cell2mat(OSI_global_perOri{g}{o}(select_fov)'));
                
                
                for tt = 1 : length(threshold)
                    if length(select_fov) < 2
                        tmp = cell2mat(OSI_perOri{g}{o}(select_fov)');
                        percent_perOri_animal{1}{tt}{g}(ff, o, 1) = sum(tmp > threshold(tt))/length(tmp);
                        percent_perOri_animal{1}{tt}{g}(ff, o, 2) = sum(tmp > threshold(tt))/(cell_no{t, g}(3, select_fov(1)));
                        tmp = cell2mat(OSI_fit_perOri{g}{o}(select_fov)');
                        percent_perOri_animal{2}{tt}{g}(ff, o, 1) = sum(tmp > threshold(tt))/length(tmp);
                        percent_perOri_animal{2}{tt}{g}(ff, o, 2) = sum(tmp > threshold(tt))/(cell_no{t, g}(3, select_fov(1)));
                        tmp = cell2mat(OSI_global_perOri{g}{o}(select_fov)');
                        percent_perOri_animal{3}{tt}{g}(ff, o, 1) = sum(tmp > threshold(tt))/length(tmp);
                        percent_perOri_animal{3}{tt}{g}(ff, o, 2) = sum(tmp > threshold(tt))/(cell_no{t, g}(3, select_fov(1)));
                    else
                        tmp = cell2mat(OSI_perOri{g}{o}(select_fov)');
                        percent_perOri_animal{1}{tt}{g}(ff, o, 1) = sum(tmp > threshold(tt))/length(tmp);
                        percent_perOri_animal{1}{tt}{g}(ff, o, 2) = sum(tmp > threshold(tt))/(cell_no{t, g}(3, select_fov(1)) + cell_no{t, g}(3, select_fov(2)));
                        tmp = cell2mat(OSI_fit_perOri{g}{o}(select_fov)');
                        percent_perOri_animal{2}{tt}{g}(ff, o, 1) = sum(tmp > threshold(tt))/length(tmp);
                        percent_perOri_animal{2}{tt}{g}(ff, o, 2) = sum(tmp > threshold(tt))/(cell_no{t, g}(3, select_fov(1)) + cell_no{t, g}(3, select_fov(2)));
                        tmp = cell2mat(OSI_global_perOri{g}{o}(select_fov)');
                        percent_perOri_animal{3}{tt}{g}(ff, o, 1) = sum(tmp > threshold(tt))/length(tmp);
                        percent_perOri_animal{3}{tt}{g}(ff, o, 2) = sum(tmp > threshold(tt))/(cell_no{t, g}(3, select_fov(1)) + cell_no{t, g}(3, select_fov(2)));
                    end
                end
                
                
                %                 count = count + 1;
            end
            
        end
    end
    
end



%% compare the percentage of cells with (global) DSI/OSI larger than thresholds across all directions (animal-based)
clear tmp
for cc = 3
    for tt = 1 : length(threshold)
        h = figure; set(h, 'position', [200 100 350 600])
        for g = 1 : length(tag)
            tmp{g} = sum(percent_perDir_animal{cc}{tt}{g}(:, :, 2), 2);
            percent_selective_animal{cc}{g}(:, tt) = tmp{g}(~isnan(tmp{g}));
            scatter(g*ones(size(percent_selective_animal{cc}{g}(:, tt))), percent_selective_animal{cc}{g}(:, tt), 'MarkerEdgeColor', colorVector(g, :), 'LineWidth', 2); hold on
        end
        [hh12(tt), p12(tt)] = ranksum(percent_selective_animal{cc}{1}(:, tt), percent_selective_animal{cc}{2}(:, tt));
        col=@(x)reshape(x,numel(x),1);
        boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)), 'whisker', 2,varargin{:});
        boxplot2(tmp, 'Width', .6, 'Colors', colorVector, 'Labels', tag, 'whisker', 1000)
        set(findobj(gca,'type','line'),'linew',4)
        set(gca,'linew',4, 'FontSize', 22)
        ylim([0 1])
        if tt == 6
            ylim([0 0.7])
        end
        box off
        title(['th=', num2str(threshold(tt)), ' p=', num2str(hh12(tt))])
        saveas(h, [save_tag, 'totalSelectivePercent_cc', num2str(cc), '_th', num2str(tt), '_animal1.png'])
    end
end


% orientation
clear tmp
for cc = 3
    for tt = 1 : length(threshold)
        h = figure; set(h, 'position', [200 100 350 600])
        for g = 1 : length(tag)
            tmp{g} = sum(percent_perOri_animal{cc}{tt}{g}(:, :, 2), 2);
            percent_selective_animal_os{cc}{g}(:, tt) = tmp{g}(~isnan(tmp{g}));
            scatter(g*ones(size(percent_selective_animal_os{cc}{g}(:, tt))), percent_selective_animal_os{cc}{g}(:, tt), 'MarkerEdgeColor', colorVector(g, :), 'LineWidth', 2); hold on
        end
        [hh13(tt), p13(tt)] = ranksum(percent_selective_animal_os{cc}{1}(:, tt), percent_selective_animal_os{cc}{2}(:, tt));
        col=@(x)reshape(x,numel(x),1);
        boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)), 'whisker', 2,varargin{:});
        boxplot2(tmp, 'Width', .6, 'Colors', colorVector, 'Labels', tag, 'whisker', 1000)
        set(findobj(gca,'type','line'),'linew',4)
        set(gca,'linew',4, 'FontSize', 22)
        ylim([0 1])
        box off
        title(['th=', num2str(threshold(tt)), ' p=', num2str(hh12(tt))])
        saveas(h, [save_tag, 'totalSelectivePercent_cc', num2str(cc), '_th', num2str(tt), '_animal1_os.png'])
    end
end



%% compare percentage of cells with (global) DSI in certain directions across mutliple thresholds (animal-based)
clear tmp
for cc = 3
    
    %     % normalized percentage
    %     figure;
    %     for g = 1 : length(tag)
    %         for tt = 1 : length(threshold)
    %             normalized_percent{g}(:, tt) = percent_perDir_animal{cc}{tt}{g}(:, 5, 2)./ sum(percent_perDir_animal{cc}{tt}{g}(:, :, 2), 2);
    %             %             scatter(tt*ones(size(normalized_percent{g}(:, tt))), normalized_percent{g}(:, tt), 'MarkerEdgeColor', colorVector(g, :)); hold on
    %             tmp = normalized_percent{g}(:, tt);
    %             mean_percent(g, tt) = mean(tmp(~isnan(tmp)));
    %             sem_percent(g, tt) = std(tmp(~isnan(tmp))) / sqrt(sum(~isnan(tmp)));
    %         end
    %         %         validId = find(isnan(normalized_percent{g}(:, 1)) == 0);
    %         %         normalized_percent{g} = normalized_percent{g}(validId, :);
    %
    % %         errorbar(threshold, mean_percent(g, :), sem_percent(g, :), 'color', colorVector(g, :)); hold on
    %     end
    %     ylim([0 0.7]); xlim([0 0.55])
    %     tmp = [];
    %     for tt = 1 : length(threshold)
    %         tmp = [tmp, normalized_percent{1}(:, tt)];
    %     end
    %     tmp2 = [];
    %     for tt = 1 : length(threshold)
    %         tmp2 = [tmp2, normalized_percent{2}(:, tt)];
    %     end
    %     [aaa bbb] = ranksum(tmp(:), tmp2(:))
    
    
    
    % percentage
    for d = 1:8
        h = figure;
        for g = 1 : length(tag)
            for tt = 1 : length(threshold)
                tmp = percent_perDir_animal{cc}{tt}{g}(:, d, 2);
                mean_percent(g, tt) = mean(tmp(~isnan(tmp)));
                sem_percent(g, tt) = std(tmp(~isnan(tmp))) / sqrt(sum(~isnan(tmp)));
                
                if tt > 1
                    start_p(1, :) = [threshold(tt-1), mean_percent(g, tt-1) + sem_percent(g, tt-1)];
                    start_p(2, :) = [threshold(tt-1), mean_percent(g, tt-1) - sem_percent(g, tt-1)];
                    end_p(1, :) = [threshold(tt), mean_percent(g, tt) + sem_percent(g, tt)];
                    end_p(2, :) = [threshold(tt), mean_percent(g, tt) - sem_percent(g, tt)];
                    fill([start_p(1, 1), end_p(1, 1), end_p(1, 1), start_p(1, 1)], [start_p(1, 2), end_p(1, 2), end_p(2, 2), start_p(2, 2)], colorVector(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on;
                    line(threshold(tt-1:tt), [mean_percent(g, tt-1), mean_percent(g, tt)], 'Color', colorVector(g, :), 'linewidth', 5); hold on
                end
            end
        end
        tmp = [];
        for tt = 1 : length(threshold)
            tmp = [tmp, percent_perDir_animal{cc}{tt}{2}(:, d, 2)];
        end
        tmp2 = [];
        for tt = 1 : length(threshold)
            tmp2 = [tmp2, percent_perDir_animal{cc}{tt}{1}(:, d, 2)];
        end
        [p10(d) hh10(d)] = ranksum(tmp(:), tmp2(:));
        title(['d', num2str(d), ' p=', num2str(p10(d))])
        box off
        set(gca, 'XTick', 0:0.1:0.5, 'linew', 4, 'FontSize', 22)
        ylim([0 0.3]);
        saveas(h, [save_tag,'normalizedProportion_allTh_d', num2str(d), '.png'])
    end
    
end


% orientation
clear tmp
for cc = 3
    % percentage
    for o = 1:4
        h = figure;
        for g = 1 : length(tag)
            for tt = 1 : length(threshold)
                tmp = percent_perOri_animal{cc}{tt}{g}(:, o, 2);
                mean_percent(g, tt) = mean(tmp(~isnan(tmp)));
                sem_percent(g, tt) = std(tmp(~isnan(tmp))) / sqrt(sum(~isnan(tmp)));
                
                if tt > 1
                    start_p(1, :) = [threshold(tt-1), mean_percent(g, tt-1) + sem_percent(g, tt-1)];
                    start_p(2, :) = [threshold(tt-1), mean_percent(g, tt-1) - sem_percent(g, tt-1)];
                    end_p(1, :) = [threshold(tt), mean_percent(g, tt) + sem_percent(g, tt)];
                    end_p(2, :) = [threshold(tt), mean_percent(g, tt) - sem_percent(g, tt)];
                    fill([start_p(1, 1), end_p(1, 1), end_p(1, 1), start_p(1, 1)], [start_p(1, 2), end_p(1, 2), end_p(2, 2), start_p(2, 2)], colorVector(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on;
                    line(threshold(tt-1:tt), [mean_percent(g, tt-1), mean_percent(g, tt)], 'Color', colorVector(g, :), 'linewidth', 5); hold on
                end
            end
        end
        tmp = [];
        for tt = 1 : length(threshold)
            tmp = [tmp, percent_perOri_animal{cc}{tt}{2}(:, o, 2)];
        end
        tmp2 = [];
        for tt = 1 : length(threshold)
            tmp2 = [tmp2, percent_perOri_animal{cc}{tt}{1}(:, o, 2)];
        end
        
        % exclude thresholds that are too high
        max_thresh1 = find(sum(tmp, 2)>0);
        max_thresh2 = find(sum(tmp2, 2)>0);
        max_thresh = max(max_thresh1(end), max_thresh2(end));
        tmp = tmp(:, 1:max_thresh);
        tmp2 = tmp2(:, 1:max_thresh);
        [pp14(o) hh14(o)] = ranksum(tmp(:), tmp2(:));
        title(['o', num2str(o), ' p=', num2str(pp14(o))])
        box off
        set(gca, 'XTick', 0:0.1:0.5, 'linew', 4, 'FontSize', 22)
        ylim([0 0.4]);
        saveas(h, [save_tag,'normalizedProportion_allTh_o', num2str(o), '.png'])
    end
end

%% compare the percentage of cells with (global) DSI/OSI larger than thresholds in certain direction (animal- and FOV-based), rose plot
cc_tag = {'DSI', 'DSI_fit', 'DSI_global'};
ds_tag = {'DS', 'DS_fit', 'DS_global'};
os_tag = {'OS', 'OS_fit', 'OS_global'};

for isanimal = [1]
    for cc = [3]
        for tt = 1 : length(threshold)
            
            h = figure;
            for g = 1 : length(tag)
                if isanimal
                    data = percent_perDir_animal{cc}{tt}{g}(:, :, 2);
                else
                    data = percent_perDir{cc}{tt}{g}(:, :, 2);
                end
                data(isnan(data)) = 0;
                proportion_perDir{cc}(g, :, tt) = mean(data);
                %             proportion_perDir{cc}(g, :, tt) = mean(data)/max(mean(data));
                
                r_vector = ((0 : 1/4 : 2) - 1/2) * pi;
                %             sem = std(data) ./ sqrt(size(data, 1)) / max(mean(data));
                sem = std(data) ./ sqrt(size(data, 1));
                for d = 1 : 8
                    start_p(1, :) = [r_vector(d), proportion_perDir{cc}(g, d, tt) + sem(d)];
                    start_p(2, :) = [r_vector(d), proportion_perDir{cc}(g, d, tt) - sem(d)];
                    if d < 8
                        end_p(1, :) = [r_vector(d+1), proportion_perDir{cc}(g, d+1, tt) + sem(d+1)];
                        end_p(2, :) = [r_vector(d+1), proportion_perDir{cc}(g, d+1, tt) - sem(d+1)];
                    else
                        end_p(1, :) = [r_vector(1), proportion_perDir{cc}(g, 1, tt) + sem(1)];
                        end_p(2, :) = [r_vector(1), proportion_perDir{cc}(g, 1, tt) - sem(1)];
                    end
                    
                    [x, y] = pol2cart([start_p(1:2, 1); end_p(1:2, 1)], [start_p(1:2, 2); end_p(1:2, 2)]);
                    fill(x([1 2 4 3]), y([1 2 4 3]), colorVector(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on;
                    avg_p = [(((0 : 1/4 : 2) - 1/2) * pi)', [proportion_perDir{cc}(g, :, tt), proportion_perDir{cc}(g, 1, tt)]'];
                    [avg_x, avg_y] = pol2cart(avg_p(:, 1), avg_p(:, 2));
                    line(avg_x, avg_y, 'Color', colorVector(g, :), 'linewidth', 5); hold on
                    hold on
                    
                    if g == 1
                        circle_r = max(proportion_perDir{cc}(g, :, tt)+sem);
                    else
                        circle_r = max(circle_r, max(proportion_perDir{cc}(g, :, tt)+sem));
                    end
                    
                    if tt == 6
                        circle_r = 0.15;
                    end
                    %                 if g == 2
                    myDrawCircle(0, 0, circle_r); hold on
                    %                 if tt == 3
                    %                     myDrawCircle(0, 0, 0.1, 1); hold on
                    %                 end
                    %                 end
                    line([-circle_r, circle_r], [0 0], 'color', [129 129 129]/255, 'lineWidth', 2, 'lineStyle', '--'); hold on
                    line([0 0], [-circle_r, circle_r], 'color', [129 129 129]/255, 'lineWidth', 2, 'lineStyle', '--'); hold on
                    set(gca,'linew',2, 'FontSize', 14, 'XColor', 'none','YColor','none');
                    title([cc_tag{cc}, ' th', num2str(tt), ' circle=', num2str(circle_r)]); axis image
                    
                end
                
            end
            saveas(h, [save_tag, cc_tag{cc}, '_percentPerDirection_th', num2str(tt), '_animal', num2str(isanimal), '.png'])
            
            
            
            % orientation
            h = figure;
            for g = 1 : length(tag)
                if isanimal
                    data = percent_perOri_animal{cc}{tt}{g}(:, :, 2);
                else
                    data = percent_perOri{cc}{tt}{g}(:, :, 2);
                end
                data(isnan(data)) = 0;
                proportion_perOri{cc}(g, :, tt) = mean(data);
                
                r_vector = ((0 : 1/4 : 2) - 1/2) * pi;
                sem = std(data) ./ sqrt(size(data, 1));
                for d = 1 : 3
                    start_p(1, :) = [d, proportion_perOri{cc}(g, d, tt) + sem(d)];
                    start_p(2, :) = [d, proportion_perOri{cc}(g, d, tt) - sem(d)];
                    
                    end_p(1, :) = [d+1, proportion_perOri{cc}(g, d+1, tt) + sem(d+1)];
                    end_p(2, :) = [d+1, proportion_perOri{cc}(g, d+1, tt) - sem(d+1)];
                    
                    fill([d d+1 d+1 d], [proportion_perOri{cc}(g, d, tt) + sem(d), proportion_perOri{cc}(g, d+1, tt) + sem(d+1), ...
                        proportion_perOri{cc}(g, d+1, tt) - sem(d+1), proportion_perOri{cc}(g, d, tt) - sem(d)], colorVector(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on;
                end
                
                line(1:4, proportion_perOri{cc}(g, :, tt), 'Color', colorVector(g, :), 'linewidth', 5); hold on
                xlim([0 5])
                if tt == 3
                    ylim([0 0.15])
                end
                %                 ylim([0 0.15])
            end
            hold off
            box off
            set(gca, 'XTick', 0:1:5, 'linew', 4, 'FontSize', 22)
            saveas(h, [save_tag, os_tag{cc}, '_percentPerOrientation_th', num2str(tt), '_animal', num2str(isanimal), '.png'])
            
            
            
            %         figure;
            %         for g = 1 : length(tag)
            %             fov_n = size(DSI_fit_perDir{g}{1});
            %             for f = 1 : fov_n(2)
            %                 proportion_th{cc}{g}(f, :) = DScell_no_thresh{cc}{tt}{g}(f, :) / sum(DScell_no_thresh{cc}{tt}{g}(f, :));
            %             end
            %
            %             proportion_th{cc}{g} = proportion_th{cc}{g}(~isnan(proportion_th{cc}{g}(:, d)), :);
            %             proportion_th_mean{cc}{g}(tt, :) = mean(proportion_th{cc}{g});
            %
            %             r_vector = ((0 : 1/4 : 2) - 1/2) * pi;
            %             sem = std(proportion_th{cc}{g}) ./ sqrt(size(proportion_th{cc}{g}, 1));
            %             for d = 1 : 8
            %                 start_p(1, :) = [r_vector(d), proportion_th_mean{cc}{g}(tt, d) + sem(d)];
            %                 start_p(2, :) = [r_vector(d), proportion_th_mean{cc}{g}(tt, d) - sem(d)];
            %                 if d < 8
            %                     end_p(1, :) = [r_vector(d+1), proportion_th_mean{cc}{g}(tt, d+1) + sem(d+1)];
            %                     end_p(2, :) = [r_vector(d+1), proportion_th_mean{cc}{g}(tt, d+1) - sem(d+1)];
            %                 else
            %                     end_p(1, :) = [r_vector(1), proportion_th_mean{cc}{g}(tt, 1) + sem(1)];
            %                     end_p(2, :) = [r_vector(1), proportion_th_mean{cc}{g}(tt, 1) - sem(1)];
            %                 end
            %                 %                     fill(p(:, 1), p(:, 2), colorVector(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.2); hold on;
            % %                 if ismember(cc, [3:5, 7:9]) && d > 3
            % %                     [x, y] = pol2cart([start_p(1:2, 1); end_p(1:2, 1)], [start_p(1:2, 2); end_p(1:2, 2)]);
            % %                     fill(x([1 2 4 3]), y([1 2 4 3]), colorVector(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.15); hold on;
            % %                 else
            %                     [x, y] = pol2cart([start_p(1:2, 1); end_p(1:2, 1)], [start_p(1:2, 2); end_p(1:2, 2)]);
            %                     fill(x([1 2 4 3]), y([1 2 4 3]), colorVector(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on;
            % %                 end
            %             end
            %             avg_p = [(((0 : 1/4 : 2) - 1/2) * pi)', [proportion_th_mean{cc}{g}(tt, :), proportion_th_mean{cc}{g}(tt, 1)]'];
            %             [avg_x, avg_y] = pol2cart(avg_p(:, 1), avg_p(:, 2));
            % %             if ismember(cc, [3:5, 7:9])
            % %                 line(avg_x(1:4), avg_y(1:4), 'Color', colorVector(g, :), 'linewidth', 5); hold on
            % %                 line(avg_x(4:end), avg_y(4:end), 'Color', colorVector(g, :), 'linewidth', 2, 'lineStyle', '--'); hold on
            % %             else
            %                 line(avg_x, avg_y, 'Color', colorVector(g, :), 'linewidth', 5); hold on
            % %             end
            %             hold on
            %
            %             if g == 1
            %                 circle_r = max(proportion_th_mean{cc}{g}(tt, :)+sem);
            %             else
            %                 circle_r = max(circle_r, max(proportion_th_mean{cc}{g}(tt, :)+sem));
            %             end
            %
            %             if g == 2
            %                 myDrawCircle(0, 0, circle_r); hold on
            %             end
            %             line([-circle_r, circle_r], [0 0], 'color', [129 129 129]/255, 'lineWidth', 2, 'lineStyle', '--'); hold on
            %             line([0 0], [-circle_r, circle_r], 'color', [129 129 129]/255, 'lineWidth', 2, 'lineStyle', '--'); hold on
            %             set(gca,'linew',2, 'FontSize', 14, 'XColor', 'none','YColor','none');
            %             title([cc_tag{cc}, ' th', num2str(tt), ' circle=', num2str(circle_r)]); axis image
            %         end
            
            
            
            for d = 1 : 8
                [p_perdir{cc}(d, tt, 1), h_perdir{cc}(d, tt, 1)] = ttest2(percent_perDir{cc}{tt}{1}(:, d, 1), percent_perDir{cc}{tt}{2}(:, d, 1));
                [p_perdir_animal{cc}(d, tt, 1), h_perdir_animal{cc}(d, tt, 1)] = ttest2(percent_perDir_animal{cc}{tt}{1}(:, d, 1), percent_perDir_animal{cc}{tt}{2}(:, d, 1));
                [p_perdir{cc}(d, tt, 2), h_perdir{cc}(d, tt, 2)] = ttest2(percent_perDir{cc}{tt}{1}(:, d, 2), percent_perDir{cc}{tt}{2}(:, d, 2));
                [p_perdir_animal{cc}(d, tt, 2), h_perdir_animal{cc}(d, tt, 2)] = ttest2(percent_perDir_animal{cc}{tt}{1}(:, d, 2), percent_perDir_animal{cc}{tt}{2}(:, d, 2));
            end
            
            
            for o = 1 : 4
                [p_perOri{cc}(o, tt, 1), h_perOri{cc}(o, tt, 1)] = ttest2(percent_perOri{cc}{tt}{1}(:, o, 1), percent_perOri{cc}{tt}{2}(:, o, 1));
                [p_perOri_animal{cc}(o, tt, 1), h_perOri_animal{cc}(o, tt, 1)] = ttest2(percent_perOri_animal{cc}{tt}{1}(:, o, 1), percent_perOri_animal{cc}{tt}{2}(:, o, 1));
                [p_perOri{cc}(o, tt, 2), h_perOri{cc}(o, tt, 2)] = ttest2(percent_perOri{cc}{tt}{1}(:, o, 2), percent_perOri{cc}{tt}{2}(:, o, 2));
                [p_perOri_animal{cc}(o, tt, 2), h_perOri_animal{cc}(o, tt, 2)] = ttest2(percent_perOri_animal{cc}{tt}{1}(:, o, 2), percent_perOri_animal{cc}{tt}{2}(:, o, 2));
            end
            
            
            %             % test whether saline animals/FOVs significantly prefer the
            %             % opposite direction
            %             % saline
            %             anova_tmp = []; anova_animal_tmp = [];
            %             group_tmp = []; group_animal_tmp = [];
            %             for d = 1 : 8
            %                 anova_tmp = [anova_tmp; percent_perDir{cc}{tt}{1}(:, d, 2)];
            %                 group_tmp = [group_tmp; d * ones(size(percent_perDir{cc}{tt}{1}(:, d, 2)))];
            %                 anova_animal_tmp = [anova_animal_tmp; percent_perDir_animal{cc}{tt}{1}(:, d, 2)];
            %                 group_animal_tmp = [group_animal_tmp; d * ones(size(percent_perDir_animal{cc}{tt}{1}(:, d, 2)))];
            %             end
            %             [bias_anova_p{cc}{tt}, bias_t{cc}{tt}, bias_stats{cc}{tt}] = anova1(anova_tmp, group_tmp);
            %             [bias_c{cc}{tt}, m] = multcompare(bias_stats{cc}{tt});
            %             [bias_animal_anova_p{cc}{tt}, bias_animal_t{cc}{tt}, bias_animal_stats{cc}{tt}] = anova1(anova_animal_tmp, group_animal_tmp, 'off');
            %             [bias_animal_c{cc}{tt}, m] = multcompare(bias_animal_stats{cc}{tt});
            %
            %             % gabazine
            %             anova_tmp = []; anova_animal_tmp = [];
            %             group_tmp = []; group_animal_tmp = [];
            %             for d = 1 : 8
            %                 anova_tmp = [anova_tmp; percent_perDir{cc}{tt}{2}(:, d, 2)];
            %                 group_tmp = [group_tmp; d * ones(size(percent_perDir{cc}{tt}{2}(:, d, 2)))];
            %                 anova_animal_tmp = [anova_animal_tmp; percent_perDir_animal{cc}{tt}{2}(:, d, 2)];
            %                 group_animal_tmp = [group_animal_tmp; d * ones(size(percent_perDir_animal{cc}{tt}{2}(:, d, 2)))];
            %             end
            %             [bias_anova_p2{cc}{tt}, bias_t2{cc}{tt}, bias_stats2{cc}{tt}] = anova1(anova_tmp, group_tmp);
            %             [bias_c2{cc}{tt}, m] = multcompare(bias_stats2{cc}{tt});
            %             [bias_animal_anova_p2{cc}{tt}, bias_animal_t2{cc}{tt}, bias_animal_stats2{cc}{tt}] = anova1(anova_animal_tmp, group_animal_tmp, 'off');
            %             [bias_animal_c2{cc}{tt}, m] = multcompare(bias_animal_stats2{cc}{tt});
            %
            %
            %
            %             % test whether saline animals/FOVs significantly prefer the
            %             % vertical orientation
            %             % saline
            %             anova_tmp = []; anova_animal_tmp = [];
            %             group_tmp = []; group_animal_tmp = [];
            %             for o = 1 : 4
            %                 anova_tmp = [anova_tmp; percent_perOri{cc}{tt}{1}(:, o, 2)];
            %                 group_tmp = [group_tmp; o * ones(size(percent_perOri{cc}{tt}{1}(:, o, 2)))];
            %                 anova_animal_tmp = [anova_animal_tmp; percent_perOri_animal{cc}{tt}{1}(:, o, 2)];
            %                 group_animal_tmp = [group_animal_tmp; o * ones(size(percent_perOri_animal{cc}{tt}{1}(:, o, 2)))];
            %             end
            %             [bias_ori_anova_p{cc}{tt}, bias_ori_t{cc}{tt}, bias_ori_stats{cc}{tt}] = anova1(anova_tmp, group_tmp);
            %             [bias_ori_c{cc}{tt}, m] = multcompare(bias_ori_stats{cc}{tt});
            %             % %             [bias_ori_animal_anova_p{cc}{tt}, bias_ori_animal_t{cc}{tt}, bias_ori_animal_stats{cc}{tt}] = anova1(anova_animal_tmp, group_animal_tmp, 'off');
            %             % %             [bias_ori_animal_c{cc}{tt}, m] = multcompare(bias_ori_animal_stats{cc}{tt});
            %
            %             % gabazine
            %             anova_tmp = []; anova_animal_tmp = [];
            %             group_tmp = []; group_animal_tmp = [];
            %             for o = 1 : 4
            %                 anova_tmp = [anova_tmp; percent_perOri{cc}{tt}{2}(:, o, 2)];
            %                 group_tmp = [group_tmp; o * ones(size(percent_perOri{cc}{tt}{2}(:, o, 2)))];
            %                 anova_animal_tmp = [anova_animal_tmp; percent_perOri_animal{cc}{tt}{2}(:, o, 2)];
            %                 group_animal_tmp = [group_animal_tmp; o * ones(size(percent_perOri_animal{cc}{tt}{2}(:, o, 2)))];
            %             end
            %             [bias_ori_anova_p2{cc}{tt}, bias_ori_t2{cc}{tt}, bias_ori_stats2{cc}{tt}] = anova1(anova_tmp, group_tmp);
            %             [bias_ori_c2{cc}{tt}, m] = multcompare(bias_ori_stats2{cc}{tt});
            %             [bias_ori_animal_anova_p2{cc}{tt}, bias_ori_animal_t2{cc}{tt}, bias_ori_animal_stats2{cc}{tt}] = anova1(anova_animal_tmp, group_animal_tmp, 'off');
            %             [bias_ori_animal_c2{cc}{tt}, m] = multcompare(bias_ori_animal_stats2{cc}{tt});
            %         end
            %
            %         data1 = percent_perDir{cc}{tt}{2}(:, 5, 2);
            %         data2 = percent_perDir{cc}{tt}{2}(:, [1:4 6:8], 2);
            %         [p_test, h_test] = ttest2(data1, data2(:))
            %
            %
            %         data1 = percent_perDir{cc}{tt}{1}(:, 5, 2);
            %         data2 = percent_perDir{cc}{tt}{1}(:, [1:4 6:8], 2);
            %         [p_test, h_test] = ttest2(data1, data2(:))
            
        end
    end
end





%% plot total indices for cells preferring different directions/orientations without thresholding: cdf and histograms
saline_label = {'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal'};
gabazine_label = {'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba'};
clear counts centers proportion

for cc = [1 3]
    ds_tag = {'DS', 'DS_fit', 'DS_global'};
    os_tag = {'OS', 'OS_fit', 'OS_global'};
    
    % cdf
    h = figure; set(h, 'position', [0 0 1000 500])
    for d = 1 : 8
        subplot(2, 4, d)
        for g = 1 : length(tag)
            data = total_index_perDir{cc}{g}{d};
            
            ch = cdfplot(data);
            ch.Color = colorVector(g, :);
            ch.LineWidth = 3;
            hold on
        end
        
        ylabel('proportion', 'FontSize', 12, 'FontWeight', 'Bold')
        [p7{1}(cc, d), hh7{1}(cc, d)] = kstest2(total_index_perDir{cc}{1}{d}, total_index_perDir{cc}{2}{d});
        [p7{2}(cc, d), hh7{2}(cc, d)] = ranksum(total_index_perDir{cc}{1}{d}, total_index_perDir{cc}{2}{d});
        title(num2str(hh7{1}(cc, d)), 'FontSize', 15)
        set(gca,'linew',2)
        xlim([0 1])
        box off; grid off
    end
    saveas(h, [save_tag, ds_tag{cc}, '_perDirection_cdf.png'])
    
    
    h = figure; set(h, 'position', [0 0 1000 250])
    for o = 1 : 4
        subplot(1, 4, o)
        for g = 1 : length(tag)
            data = total_index_perOri{cc}{g}{o};
            
            ch = cdfplot(data);
            ch.Color = colorVector(g, :);
            ch.LineWidth = 3;
            hold on
        end
        
        ylabel('proportion', 'FontSize', 12, 'FontWeight', 'Bold')
        [p7{1}(cc+3, o), hh7{1}(cc+3, o)] = kstest2(total_index_perOri{cc}{1}{o}, total_index_perOri{cc}{2}{o});
%         [p15{2}(cc+3, o), hh15{2}(cc+3, o)] = ttest2(total_index_perOri{cc}{1}{o}, total_index_perOri{cc}{2}{o});
        title(num2str(hh7{1}(cc+3, o)), 'FontSize', 15)
        set(gca,'linew',2)
        xlim([0 1])
        box off; grid off
    end
    saveas(h, [save_tag, os_tag{cc}, '_perOrientation_cdf.png'])
    
    
    
    
    % histograms
    anova_tmp = [];
    clear group_tmp
    group_tmp{1} = [];
    group_tmp{2} = [];
    clear tmp
    h = figure; set(h, 'position', [0 0 1000 500])
    for d = 1 : 8
        subplot(2, 4, d)
        for g = 1 : length(tag)
            data = total_index_perDir{cc}{g}{d};
            
            [counts(g, :), centers(g, :)] = hist(data, 0.05:0.1:1);
            proportion(g, :) = counts(g, :) / length(data);
            if g == 1
                facevalue = 0.3;
            else
                facevalue = 0.2;
            end
            bar(centers(g, :), proportion(g, :), 'barWidth', 0.9, 'FaceColor', colorVector(g, :),...
                'FaceAlpha', facevalue, 'EdgeColor', colorVector_line(g, :), 'LineWidth', 1.5); hold on
            tmp{(d-1)*2+g} = data;
        end
        anova_tmp = [anova_tmp; tmp{(d-1)*2+1}; tmp{(d-1)*2+2}];
        group_tmp{1} = [group_tmp{1}; d * ones(length(tmp{(d-1)*2+1}), 1); d * ones(length(tmp{(d-1)*2+2}), 1)];
        group_tmp{2} = [group_tmp{2}; repmat(saline_label(1), length(tmp{(d-1)*2+1}), 1); repmat(gabazine_label(1), length(tmp{(d-1)*2+2}), 1)];
        
        ylabel('proportion', 'FontSize', 12, 'FontWeight', 'Bold')
        title(num2str(hh7{2}(cc, d)), 'FontSize', 15)
        set(gca,'linew',2)
        xlim([0 1])
        box off; grid off
    end
    saveas(h, [save_tag, ds_tag{cc}, '_perDirection_hist.png'])
    [direction_anova_p{1}{cc}, direction_t{1}{cc}, direction_stats{1}{cc}] = anovan(anova_tmp, {group_tmp{1}, group_tmp{2}}, 'interaction');
    [direction_c{1}{cc}, m] = multcompare(direction_stats{1}{cc}, 'dimension', [1, 2], 'CType','bonferroni');
    
    
    anova_tmp = [];
    group_tmp{1} = [];
    group_tmp{2} = [];
    clear tmp
    h = figure; set(h, 'position', [0 0 1000 250])
    for o = 1 : 4
        subplot(1, 4, o)
        for g = 1 : length(tag)
            data = total_index_perOri{cc}{g}{o};
            
            [counts(g, :), centers(g, :)] = hist(data, 0.05:0.1:1);
            proportion(g, :) = counts(g, :) / length(data);
            if g == 1
                facevalue = 0.3;
            else
                facevalue = 0.2;
            end
            bar(centers(g, :), proportion(g, :), 'barWidth', 0.9, 'FaceColor', colorVector(g, :),...
                'FaceAlpha', facevalue, 'EdgeColor', colorVector_line(g, :), 'LineWidth', 1.5); hold on
            tmp{(o-1)*2+g} = data;
        end
        anova_tmp = [anova_tmp; tmp{(o-1)*2+1}; tmp{(o-1)*2+2}];
        group_tmp{1} = [group_tmp{1}; o * ones(length(tmp{(o-1)*2+1}), 1); o * ones(length(tmp{(o-1)*2+2}), 1)];
        group_tmp{2} = [group_tmp{2}; repmat(saline_label(1), length(tmp{(o-1)*2+1}), 1); repmat(gabazine_label(1), length(tmp{(o-1)*2+2}), 1)];
        ylabel('proportion', 'FontSize', 12, 'FontWeight', 'Bold')
        title(num2str(hh7{2}(cc+3, o)), 'FontSize', 15)
        set(gca,'linew',2)
        xlim([0 1])
        box off; grid off
    end
    saveas(h, [save_tag, os_tag{cc}, '_perOrientation_hist.png'])
    [orientation_anova_p{1}{cc}, t, orientation_stats{1}{cc}] = anovan(anova_tmp, {group_tmp{1}, group_tmp{2}}, 'interaction');
    [orientation_c{1}{cc}, m] = multcompare(orientation_stats{1}{cc}, 'dimension', [1, 2]);
    
end




%% plot total indices for cells preferring horizontal directions vs other directons without thresholding: cdf and histograms
saline_label = {'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal'};
gabazine_label = {'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba'};

clear counts centers proportion

for cc = 3
    ds_tag2 = {'DS_horizontal', 'DS_others'};
    
    % cdf
    for group_d = 1 : 2
        h = figure; %set(h, 'position', [0 0 1000 500])
        for g = 1 : length(tag)
            if group_d == 1
                data = [total_index_perDir{cc}{g}{5}];
            else
                data = [total_index_perDir{cc}{g}{2}; total_index_perDir{cc}{g}{3}; total_index_perDir{cc}{g}{4}; total_index_perDir{cc}{g}{6}; ...
                    total_index_perDir{cc}{g}{7}; total_index_perDir{cc}{g}{8}; total_index_perDir{cc}{g}{1}];
            end
            ch = cdfplot(data);
            ch.Color = colorVector(g, :);
            ch.LineWidth = 3;
            hold on
            if g == 1
                pre_data = data;
            end
        end
        ylabel('proportion', 'FontSize', 12, 'FontWeight', 'Bold')
        [p9{1}(cc, group_d), hh9{1}(cc, group_d)] = kstest2(pre_data, data);
        [hh9{2}(cc, group_d), p9{2}(cc, group_d)] = ranksum(pre_data, data);  % changed to non-parametric test on 11/28/20
        %         [p9{2}(cc, group_d), hh9{2}(cc, group_d)] = ttest2(pre_data, data);
        title([ds_tag{cc}, ' ', ds_tag2{group_d}, ' ', num2str(hh9{1}(cc, group_d))], 'FontSize', 15)
        set(gca,'linew',3, 'FontSize', 22); axis image
        xlim([0 1])
        box off; grid off
        saveas(h, [save_tag, ds_tag{cc}, '_', ds_tag2{group_d}, '_perDirection_cdf.png'])
        
        
        
        h = figure; %set(h, 'position', [0 0 1000 500])
        for g = 1 : length(tag)
            if group_d == 1
                data = [total_index_perDir{cc}{g}{5}];
            else
                data = [total_index_perDir{cc}{g}{2}; total_index_perDir{cc}{g}{3}; total_index_perDir{cc}{g}{4}; total_index_perDir{cc}{g}{6}; ...
                    total_index_perDir{cc}{g}{7}; total_index_perDir{cc}{g}{8};total_index_perDir{cc}{g}{1}];
            end
            median_data = median(data);
            [counts(g, :), centers(g, :)] = hist(data, 0.025:0.05:1);
            proportion(g, :) = counts(g, :) / length(data);
            if g == 1
                facevalue = 0.3;
            else
                facevalue = 0.2;
            end
            bar(centers(g, :), proportion(g, :), 'barWidth', 0.9, 'FaceColor', colorVector(g, :),...
                'FaceAlpha', facevalue, 'EdgeColor', colorVector_line(g, :), 'LineWidth', 1.5); hold on
            line(median_data * [1 1], [0 0.4], 'lineStyle', '- -', 'color', colorVector(g, :), 'lineWidth', 3); hold on % for frmd7
            %             line(median_data * [1 1], [0 0.25], 'lineStyle', '- -', 'color', colorVector(g, :), 'lineWidth', 3); hold on
        end
        ylabel('proportion', 'FontSize', 12, 'FontWeight', 'Bold')
        title([ds_tag{cc}, ' ', ds_tag2{group_d}, ' ', num2str(hh9{2}(cc, group_d))], 'FontSize', 15)
        set(gca,'linew',3, 'FontSize', 22)
        ylim([0 0.4]) % for frmd7
%         ylim([0 0.25])
        xlim([0 0.8])
        box off; grid off
        saveas(h, [save_tag, ds_tag{cc}, '_', ds_tag2{group_d}, '_perDirection_hist.png'])
    end
    
end





%% test for individual median DSI or OSI per direction/orientation
clear group_tmp

group = {'1', '1', ' ', '2', '2', ' ', '3', '3', ' ', '4', '4', ' ', '5', '5', ' ', '6', '6', ' ', '7', '7', ' ', '8', '8', ' '};
saline_label = {'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal', 'sal'};
gabazine_label = {'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba', 'gaba'};


for isanimal = 1
    for cc = 1
        ds_tag = {'DS', 'DS_fit', 'DS_global'};
        os_tag = {'OS', 'OS_fit', 'OS_global'};
        anova_tmp = [];
        group_tmp{1} = [];
        group_tmp{2} = [];
        clear tmp
        clear X
        
        % scatter
        h = figure; set(h, 'position', [0 0 1000 500])
        for d = 1 : 8
            for g = 1 : length(tag)
                if isanimal
                    data = median_index_perDir_animal{cc}{g}(:, d);
                else
                    data = median_index_perDir{cc}{g}(:, d);
                end
                scatter(((d-1)*3+g) * ones(size(data)), data, 'MarkerEdgeColor', colorVector(g, :), 'LineWidth', 2); hold on
                tmp{(d-1)*3+g} = data(~isnan(data));
            end
            anova_tmp = [anova_tmp; tmp{(d-1)*3+1}; tmp{(d-1)*3+2}];
            group_tmp{1} = [group_tmp{1}; d * ones(length(tmp{(d-1)*3+1}), 1); d * ones(length(tmp{(d-1)*3+2}), 1)];
            group_tmp{2} = [group_tmp{2}; saline_label(1:length(tmp{(d-1)*3+1}))'; gabazine_label(1:length(tmp{(d-1)*3+2}))'];
            tmp{d*3} = -1;
            ylabel('proportion', 'FontSize', 12, 'FontWeight', 'Bold')
            [p8{isanimal+1}(cc, d), hh8{isanimal+1}(cc, d)] = ranksum(tmp{(d-1)*3+1}, tmp{(d-1)*3+2});
            %             [p8{isanimal+1}(cc, d), hh8{isanimal+1}(cc, d)] = ttest2(tmp{(d-1)*3+1}, tmp{(d-1)*3+2});
            title(ds_tag{cc})
            %             if d == 1 && g == 1
            %                 X(:, 1) = data;
            %                 X(:, 2) = g * ones(size(data));
            %                 X(:, 3) = d * ones(size(data));
            %                 X(:, 4) = repmat([1:27]', 8, 1);
            %             else
            %                 X(:, 1) = anova_tmp;
            %                 X(:, 2) = repmat([ones(13, 1); 2*ones(14, 1)], 8, 1);
            %                 X(:, 3) = group_tmp{1};
            %                 X(:, 4) = repmat([1:27]', 8, 1);
            %             end
        end
        col=@(x)reshape(x,numel(x),1);
        boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)), 'whisker', 2,varargin{:});
        boxplot2(tmp, 'Width', .6, 'Colors', colorVector([1 2 1 1 2 1 1 2 1 1 2 1 1 2 1 1 2 1 1 2 1 1 2 1], :), 'Labels', group, 'whisker', 1000)
        set(findobj(gca,'type','line'),'linew',4)
        set(gca,'linew',4, 'FontSize', 22)
        ylim([0 0.9])
        box off
        saveas(h, [save_tag, ds_tag{cc}, '_perDirection_median_animal', num2str(isanimal), '.png'])
        
        if isanimal
            % non-parametric 11/28/20
            %             [direction_anova_p_animal{2}{cc}, direction_t_animal{2}{cc}, direction_stats_animal{2}{cc}] = friedman(anova_tmp, {group_tmp{1}, group_tmp{2}}, 'interaction');
            %             [direction_c_animal{2}{cc}, m] = multcompare(direction_stats_animal{2}{cc}, 'dimension', [1, 2]);
            % %             [direction_anova_p_animal{2}{cc}, direction_t_animal{2}{cc}, direction_stats_animal{2}{cc}] = anovan(anova_tmp, {group_tmp{1}, group_tmp{2}}, 'interaction');
            %             [direction_c_animal{2}{cc}, m] = multcompare(direction_stats_animal{2}{cc}, 'dimension', [1, 2]);
        else
            [direction_anova_p{2}{cc}, direction_t{2}{cc}, direction_stats{2}{cc}] = anovan(anova_tmp, {group_tmp{1}, group_tmp{2}}, 'interaction');
            [direction_c{2}{cc}, m] = multcompare(direction_stats{2}{cc}, 'dimension', [1, 2]);
        end
        %         [SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(X, 0);
        %         SSQs
        
        
        
        
        
        
        h = figure; set(h, 'position', [0 0 1000 500])
        anova_tmp = [];
        group_tmp{1} = [];
        group_tmp{2} = [];
        clear tmp
        for o = 1 : 4
            for g = 1 : length(tag)
                if isanimal
                    data = median_index_perOri_animal{cc}{g}(:, o);
                else
                    data = median_index_perOri{cc}{g}(:, o);
                end
                scatter(((o-1)*3+g) * ones(size(data)), data, 'MarkerEdgeColor', colorVector(g, :), 'LineWidth', 2); hold on
                tmp{(o-1)*3+g} = data;
            end
            
            anova_tmp = [anova_tmp; tmp{(o-1)*3+1}; tmp{(o-1)*3+2}];
            group_tmp{1} = [group_tmp{1}; o * ones(length(tmp{(o-1)*3+1}), 1); o * ones(length(tmp{(o-1)*3+2}), 1)];
            group_tmp{2} = [group_tmp{2}; saline_label(1:length(tmp{(o-1)*3+1}))'; gabazine_label(1:length(tmp{(o-1)*3+2}))'];
            tmp{o*3} = -1;
            
            ylabel('proportion', 'FontSize', 12, 'FontWeight', 'Bold')
            [p8{isanimal+1}(cc+3, o), hh8{isanimal+1}(cc+3, o)] = ttest2(tmp{(o-1)*3+1}, tmp{(o-1)*3+2});
            title(os_tag{cc})
        end
        col=@(x)reshape(x,numel(x),1);
        boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)), 'whisker', 2,varargin{:});
        boxplot2(tmp, 'Width', .6, 'Colors', colorVector([1 2 1 1 2 1 1 2 1 1 2 1], :), 'Labels', group(1:12), 'whisker', 1000)
        set(findobj(gca,'type','line'),'linew',4)
        set(gca,'linew',4, 'FontSize', 22)
        ylim([0 0.6])
        box off
        saveas(h, [save_tag, os_tag{cc}, '_perOrientation_median_animal', num2str(isanimal), '.png'])
        
        if isanimal
            [orientation_anova_p_animal{2}{cc}, orientation_t_animal, orientation_stats_animal{2}{cc}] = anovan(anova_tmp, {group_tmp{1}, group_tmp{2}}, 'interaction');
            [orientation_c_animal{2}{cc}, m] = multcompare(orientation_stats_animal{2}{cc}, 'dimension', [1, 2]);
        else
            [orientation_anova_p{2}{cc}, orientation_t, orientation_stats{2}{cc}] = anovan(anova_tmp, {group_tmp{1}, group_tmp{2}}, 'interaction');
            [orientation_c{2}{cc}, m] = multcompare(orientation_stats{2}{cc}, 'dimension', [1, 2]); %'CType','bonferroni'
        end
        clear X
        %         X(:, 1) = anova_tmp;
        %         X(:, 2) = repmat([ones(13, 1); 2*ones(14, 1)], 4, 1);
        %         X(:, 3) = group_tmp{1};
        %         X(:, 4) = repmat([1:27]', 4, 1);
        %         [SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(X, 0);
        %         SSQs
        
        
        % plot gDSI for all directions
        h = figure;
        for g = 1 : length(tag)
            for d = 1 : 8
                if isanimal
                    data = median_index_perDir_animal{cc}{g}(:, d);
                else
                    data = median_index_perDir{cc}{g}(:, d);
                end
                data = data(~isnan(data));
                mean_medianIndex{cc}(g, d, tt) = mean(data);
                
                r_vector = ((0 : 1/4 : 2) - 1/2) * pi;
                sem(d) = std(data) ./ sqrt(size(data, 1));
            end
            for d = 1 : 8
                start_p(1, :) = [r_vector(d), mean_medianIndex{cc}(g, d, tt) + sem(d)];
                start_p(2, :) = [r_vector(d), mean_medianIndex{cc}(g, d, tt) - sem(d)];
                if d < 8
                    end_p(1, :) = [r_vector(d+1), mean_medianIndex{cc}(g, d+1, tt) + sem(d+1)];
                    end_p(2, :) = [r_vector(d+1), mean_medianIndex{cc}(g, d+1, tt) - sem(d+1)];
                else
                    end_p(1, :) = [r_vector(1), mean_medianIndex{cc}(g, 1, tt) + sem(1)];
                    end_p(2, :) = [r_vector(1), mean_medianIndex{cc}(g, 1, tt) - sem(1)];
                end
                
                [x, y] = pol2cart([start_p(1:2, 1); end_p(1:2, 1)], [start_p(1:2, 2); end_p(1:2, 2)]);
                fill(x([1 2 4 3]), y([1 2 4 3]), colorVector(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on;
                avg_p = [(((0 : 1/4 : 2) - 1/2) * pi)', [mean_medianIndex{cc}(g, :, tt), mean_medianIndex{cc}(g, 1, tt)]'];
                [avg_x, avg_y] = pol2cart(avg_p(:, 1), avg_p(:, 2));
                line(avg_x, avg_y, 'Color', colorVector(g, :), 'linewidth', 5); hold on
                hold on
                
                if g == 1
                    circle_r = max(mean_medianIndex{cc}(g, :, tt)+sem);
                else
                    circle_r = max(circle_r, max(mean_medianIndex{cc}(g, :, tt)+sem));
                end
                
                if cc == 3
                    circle_r = 0.4;
                elseif cc == 1
                    circle_r = 0.65;
                end
                
                if g == 2
                    myDrawCircle(0, 0, circle_r); hold on
                end
                line([-circle_r, circle_r], [0 0], 'color', [129 129 129]/255, 'lineWidth', 2, 'lineStyle', '--'); hold on
                line([0 0], [-circle_r, circle_r], 'color', [129 129 129]/255, 'lineWidth', 2, 'lineStyle', '--'); hold on
                set(gca,'linew',2, 'FontSize', 14, 'XColor', 'none','YColor','none');
                title([ds_tag{cc}, ' th', num2str(tt), ' circle=', num2str(circle_r)]); axis image
                
            end
            
        end
        saveas(h, [save_tag, ds_tag{cc}, '_median_th', num2str(tt), '_animal', num2str(isanimal),'.png'])
        
        
        
        
        % plot gOSI for all directions
        h = figure;
        for g = 1 : length(tag)
            for o = 1 : 4
                data = median_index_perOri{cc}{g}(:, o);
                data = data(~isnan(data));
                mean_medianIndex_ori{cc}(g, o, tt) = mean(data);
                
                r_vector = ((0 : 1/4 : 2) - 1/2) * pi;
                sem(d) = std(data) ./ sqrt(size(data, 1));
            end
            for d = 1 : 3
                start_p(1, :) = [d, mean_medianIndex_ori{cc}(g, d, tt) + sem(d)];
                start_p(2, :) = [d, mean_medianIndex_ori{cc}(g, d, tt) - sem(d)];
                
                end_p(1, :) = [d+1, mean_medianIndex_ori{cc}(g, d+1, tt) + sem(d+1)];
                end_p(2, :) = [d+1, mean_medianIndex_ori{cc}(g, d+1, tt) - sem(d+1)];
                
                fill([d d+1 d+1 d], [mean_medianIndex_ori{cc}(g, d, tt) + sem(d), mean_medianIndex_ori{cc}(g, d+1, tt) + sem(d+1), ...
                    mean_medianIndex_ori{cc}(g, d+1, tt) - sem(d+1), mean_medianIndex_ori{cc}(g, d, tt) - sem(d)], colorVector(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on;
            end
            
            line(1:4, mean_medianIndex_ori{cc}(g, :, tt), 'Color', colorVector(g, :), 'linewidth', 5); hold on
            xlim([0 5])
            ylim([0 0.2])
            if cc == 1
                ylim([0 0.4])
            end
        end
        hold off
        box off
        set(gca, 'XTick', 0:1:5, 'linew', 4, 'FontSize', 22)
        saveas(h, [save_tag, os_tag{cc}, '_median_th', num2str(tt), '.png'])
        
    end
end



% saline
anova_tmp = [];
group_tmp = [];
for d = 1 : 8
    anova_tmp = [anova_tmp; median_index_perDir{cc}{1}(:, d)];
    group_tmp = [group_tmp; d * ones(size(median_index_perDir{cc}{1}(:, d)))];
end
[bias_median_anova_p{cc}, bias_median_t{cc}, bias_median_stats{cc}] = anova1(anova_tmp, group_tmp);
[bias_median_c{cc}, m] = multcompare(bias_median_stats{cc});


% gabazine
anova_tmp = [];
group_tmp = [];
for d = 1 : 8
    anova_tmp = [anova_tmp; median_index_perDir{cc}{2}(:, d)];
    group_tmp = [group_tmp; d * ones(size(median_index_perDir{cc}{2}(:, d)))];
end
[bias_median_anova_p2{cc}, bias_median_t2{cc}, bias_median_stats2{cc}] = anova1(anova_tmp, group_tmp);
[bias_median_c2{cc}, m] = multcompare(bias_median_stats2{cc});
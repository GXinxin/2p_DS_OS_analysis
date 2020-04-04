%% responding cell percentage

for cc = 1 : 5
    h = figure;
    set(h, 'position', [0 0 1500 450])
    
    for t = 1 : 3
        subplot(1, 3, t)
        
        for g = 1:length(tag)
            switch cc
                case 1
                    data = resp_Percent_all{t, g};
                    cc_tag = 'responding';
                case 2
                    data = percent_DS{t, g};
                    cc_tag = 'DSI';
                case 3
                    data = percent_DS_fit{t, g};
                    cc_tag = 'DSI_fit';
                case 4
                    data = percent_OS{t, g};
                    cc_tag = 'OSI';
                case 5
                    data = percent_OS_fit{t, g};
                    cc_tag = 'OSI_fit';
%                 case 6
%                     data = percent_OS_DS{t, g};
%                     cc_tag = 'OS_DS';
%                 case 7
%                     data = percent_OS_only{t, g};
%                     cc_tag = 'OS_only';
%                 case 8
%                     data = percent_OS_DS_fit{t, g};
%                     cc_tag = 'OS_DS_fit';
%                 case 9
%                     data = percent_OS_only_fit{t, g};
%                     cc_tag = 'OS_only_fit';
            end
            
            if g == 1
                pre_data = data;
            end
            
            if cc == 1
                scatter(g * ones(1, length(data)), data, 'MarkerEdgeColor', colorVector(g, :)', 'LineWidth', 2)
                hold on
            else
                scatter(g * ones(1, length(data(1, :))), data(1, :), 'MarkerEdgeColor', colorVector(g, :)', 'LineWidth', 2)
                hold on
                scatter((g+2) * ones(1, length(data(2, :))), data(2, :), 'MarkerEdgeColor', colorVector(g, :)', 'LineWidth', 2)
                hold on
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
        boxplot2(tmp, 'Width', .6, 'Colors', colors, 'Labels', ftag, 'whisker', 1000)
        
        title([cc_tag, ' t', num2str(t)])
        ylabel('proportion')
        set(findobj(gca,'type','line'),'linew',3)
        set(gca,'linew',4, 'FontSize', 14)
        if cc == 1
            xlim([0 length(tag)+1])
        else
            xlim([0 2*length(tag)+1])
        end
        ylim([-0.1 1])
        box off
        saveas(h, [save_tag, cc_tag, '_percent.png'])
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





%% preferred direction/orientation distribution

% compute proportion for fitted data
figure;
for th = 1 : 2
    for t = 1 : 3
        for g = 1:length(tag)
            h = polarhistogram(deg2rad(DS_direction_fit_total{t, th, g}), 'BinEdges', -1/8*pi : 1/4*pi : 15/8*pi);
            DS_fit_counts{th, t}(g, :) = h.BinCounts;
            DS_fit_proportion{th, t}(g, :) = DS_fit_counts{th, t}(g, :) / length(DS_direction_fit_total{t, th, g});
            
            data = OS_direction_fit_total{t, th, g};
            data = [data, data + 180];
            h = polarhistogram(deg2rad(data), 'BinEdges', -1/8*pi : 1/4*pi : 15/8*pi);
            OS_fit_counts{th, t}(g, :) = h.BinCounts;
            OS_fit_proportion{th, t}(g, :) = 2 * OS_fit_counts{th, t}(g, :) / length(data);
            
            
%             data = fit_OS_DS_direction_total{t, th, g};
%             data = [data, data + 180];
%             h = polarhistogram(deg2rad(data), 'BinEdges', -1/8*pi : 1/4*pi : 15/8*pi);
%             fit_OS_DS_counts{th, t}(g, :) = h.BinCounts;
%             fit_OS_DS_proportion{th, t}(g, :) = 2 * fit_OS_DS_counts{th, t}(g, :) / length(data);
%             
%             data = fit_OS_only_direction_total{t, th, g};
%             data = [data, data + 180];
%             h = polarhistogram(deg2rad(data), 'BinEdges', -1/8*pi : 1/4*pi : 15/8*pi);
%             fit_OS_only_counts{th, t}(g, :) = h.BinCounts;
%             fit_OS_only_proportion{th, t}(g, :) = 2 * fit_OS_only_counts{th, t}(g, :) / length(data);
%             
%             data = fit_OS_DS_direction_direction_total{t, th, g};
%             h = polarhistogram(deg2rad(data), 'BinEdges', -1/8*pi : 1/4*pi : 15/8*pi);
%             fit_OS_DS_direction_counts{th, t}(g, :) = h.BinCounts;
%             fit_OS_DS_direction_proportion{th, t}(g, :) = fit_OS_DS_direction_counts{th, t}(g, :) / length(data);
        end
    end
end


for cc = 1 : 4
    count = 1;
    h = figure;
    set(h, 'position', [0 0 1000 1000])
    
    for th = 1 : 2
        
        for t = 1 : 3
            subplot(2, 3, count)
            
            for g = 1:length(tag)
                switch cc
                    case 1
                        data = DS_direction_total{t, th, g};
                        interval = 1:8;
                        cc_tag = 'DS';
                    case 2
                        proportion_d(g, :) = DS_fit_proportion{th, t}(g, :);
                        cc_tag = 'DS_fit';
                    case 3
                        data = OS_direction_total{t, th, g};
                        data = [data; data + 4];
                        interval = 1:8;
                        cc_tag = 'OS';
                    case 4
                        proportion_d(g, :) = OS_fit_proportion{th, t}(g, :);
                        cc_tag = 'OS_fit';
                end
                
                if cc == 1 || cc == 3
                    [counts_d(g, :), centers_d(g, :)] = hist(data, interval);
                    proportion_d(g, :) = counts_d(g, :) / length(data);
                    if cc == 3
                        proportion_d(g, :) = proportion_d(g, :) * 2;
                    end
                end
                polarplot(((0 : 1/4 : 2) - 1/2) * pi, [proportion_d(g, :), proportion_d(g, 1)], 'Color', colorVector(g, :), 'linewidth', 4)
                hold on
            end
            
            title([cc_tag, ' th', num2str(th), ' t', num2str(t)])
            set(gca,'linew',2, 'FontSize', 14)
            count = count + 1;
        end
    end
    
    saveas(h, [save_tag, cc_tag, '_rose.png'])
end




%% compute proportions of cells for each direction (FOV-based)
for g = 1 : 2
    for ii = 1 : length(DS_direction_individual{2, 1, g})
        % DS
        for d = 1 : 8
            direction_count_n{1}(ii, d, g) = sum(DS_direction_individual{2, 1, g}{ii} == d) / length(DS_direction_individual{2, 1, g}{ii});
        end
        
        % DS fit
        data = DS_direction_individual_fit{2, 1, g}{ii};
        tmp_h = polarhistogram(deg2rad(data), 'BinEdges', -1/8*pi : 1/4*pi : 15/8*pi);
        direction_count_n{2}(ii, :, g) = tmp_h.BinCounts / length(data);
        
        % OS
        data = [OS_direction_individual{2, 1, g}{ii}; OS_direction_individual{2, 1, g}{ii} + 4];
        for d = 1 : 8
            direction_count_n{3}(ii, d, g) = sum(data == d) / length(data);
        end
        
        % OS fit
        data = [OS_direction_individual_fit{2, 1, g}{ii}; OS_direction_individual_fit{2, 1, g}{ii} + 180];
        tmp_h = polarhistogram(deg2rad(data), 'BinEdges', -1/8*pi : 1/4*pi : 15/8*pi);
        direction_count_n{4}(ii, :, g) = tmp_h.BinCounts / length(data);
        
%         % OSDS
%         data = [OS_DS_direction_individual{2, 1, g}{ii}; OS_DS_direction_individual{2, 1, g}{ii} + 4];
%         for d = 1 : 8
%             direction_count_n{5}(ii, d, g) = sum(data == d) / length(data);
%         end
%         
%         % OSDS_direction
%         data = OS_DS_direction_direction_total_individual{2, 1, g}{ii};
%         for d = 1 : 8
%             direction_count_n{6}(ii, d, g) = sum(data == d) / length(data);
%         end
%         
%         % OS_only
%         data = [OS_only_direction_total_individual{2, 1, g}{ii}; OS_only_direction_total_individual{2, 1, g}{ii} + 4];
%         for d = 1 : 8
%             direction_count_n{7}(ii, d, g) = sum(data == d) / length(data);
%         end
%         
%         % OSDS fit
%         data = [OS_DS_direction_total_individual_fit{2, 1, g}{ii}; OS_DS_direction_total_individual_fit{2, 1, g}{ii} + 180];
%         tmp_h = polarhistogram(deg2rad(data), 'BinEdges', -1/8*pi : 1/4*pi : 15/8*pi);
%         direction_count_n{8}(ii, :, g) = tmp_h.BinCounts / length(data);
%         
%         % OS_only fit
%         data = [OS_only_direction_total_individual_fit{2, 1, g}{ii}; OS_only_direction_total_individual_fit{2, 1, g}{ii} + 180];
%         tmp_h = polarhistogram(deg2rad(data), 'BinEdges', -1/8*pi : 1/4*pi : 15/8*pi);
%         direction_count_n{9}(ii, :, g) = tmp_h.BinCounts / length(data);
        
%         % OS fit
%         data = OS_DS_direction_direction_total_individual_fit{2, 1, g}{ii};
%         tmp_h = polarhistogram(deg2rad(data), 'BinEdges', -1/8*pi : 1/4*pi : 15/8*pi);
%         direction_count_n{10}(ii, :, g) = tmp_h.BinCounts / length(data);
%         
        
    end
end




%% compute proportions of cells for each direction (animal-based)
for g = 1 : 2
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
                
            % DS 
            if length(select_fov) == 2
                data = [DS_direction_individual{2, 1, g}{select_fov(1)}; DS_direction_individual{2, 1, g}{select_fov(2)}];
            elseif length(select_fov) == 1
                data = DS_direction_individual{2, 1, g}{select_fov(1)};
            end
            for d = 1 : 8
                direction_count_n_animal{1}(a, d, g) = sum(data == d) / length(data);
            end


            % DS fit
            if length(select_fov) == 2
                data = [DS_direction_individual_fit{2, 1, g}{select_fov(1)}, DS_direction_individual_fit{2, 1, g}{select_fov(2)}];
            elseif length(select_fov) == 1
                data = DS_direction_individual_fit{2, 1, g}{select_fov};
            end
            tmp_h = polarhistogram(deg2rad(data), 'BinEdges', -1/8*pi : 1/4*pi : 15/8*pi);
            direction_count_n_animal{2}(a, :, g) = tmp_h.BinCounts / length(data);


            % OS
            if length(select_fov) == 2
                data = [OS_direction_individual{2, 1, g}{select_fov(1)}; OS_direction_individual{2, 1, g}{select_fov(1)} + 4; ...
                    OS_direction_individual{2, 1, g}{select_fov(2)}; OS_direction_individual{2, 1, g}{select_fov(2)} + 4];
            elseif length(select_fov) == 1
                data = [OS_direction_individual{2, 1, g}{select_fov}; OS_direction_individual{2, 1, g}{select_fov} + 4];
            end
            for d = 1 : 8
                direction_count_n_animal{3}(a, d, g) = sum(data == d) / length(data);
            end


            % OS fit
            if length(select_fov) == 2
                data = [OS_direction_individual_fit{2, 1, g}{select_fov(1)}, OS_direction_individual_fit{2, 1, g}{select_fov(1)} + 180, ...
                    OS_direction_individual_fit{2, 1, g}{select_fov(2)}, OS_direction_individual_fit{2, 1, g}{select_fov(2)} + 180];
            elseif length(select_fov) == 1
                data = [OS_direction_individual_fit{2, 1, g}{select_fov}, OS_direction_individual_fit{2, 1, g}{select_fov} + 180];
            end
            tmp_h = polarhistogram(deg2rad(data), 'BinEdges', -1/8*pi : 1/4*pi : 15/8*pi);
            direction_count_n_animal{4}(a, :, g) = tmp_h.BinCounts / length(data);


%             % OSDS
%             if length(select_fov) == 2
%                 data = [OS_DS_direction_individual{2, 1, g}{select_fov(1)}; OS_DS_direction_individual{2, 1, g}{select_fov(1)} + 4; ...
%                     OS_DS_direction_individual{2, 1, g}{select_fov(2)}; OS_DS_direction_individual{2, 1, g}{select_fov(2)} + 4];
%             elseif length(select_fov) == 1
%                 data = [OS_DS_direction_individual{2, 1, g}{select_fov}; OS_DS_direction_individual{2, 1, g}{select_fov} + 4];
%             end
%             for d = 1 : 8
%                 direction_count_n_animal{5}(a, d, g) = sum(data == d) / length(data);
%             end
% 
% 
%             % OSDS_direction
%             if length(select_fov) == 2
%                 data = [OS_DS_direction_direction_total_individual{2, 1, g}{select_fov(1)}; OS_DS_direction_direction_total_individual{2, 1, g}{select_fov(2)}];
%             elseif length(select_fov) == 1
%                 data = OS_DS_direction_direction_total_individual{2, 1, g}{select_fov(1)};
%             end
%             for d = 1 : 8
%                 direction_count_n_animal{6}(a, d, g) = sum(data == d) / length(data);
%             end
% 
% 
%             % OS_only
%             if length(select_fov) == 2
%                 data = [OS_only_direction_total_individual{2, 1, g}{select_fov(1)}; OS_only_direction_total_individual{2, 1, g}{select_fov(1)} + 4; ...
%                     OS_only_direction_total_individual{2, 1, g}{select_fov(2)}; OS_only_direction_total_individual{2, 1, g}{select_fov(2)} + 4];
%             elseif length(select_fov) == 1
%                 data = [OS_only_direction_total_individual{2, 1, g}{select_fov}; OS_only_direction_total_individual{2, 1, g}{select_fov} + 4];
%             end
%             for d = 1 : 8
%                 direction_count_n_animal{7}(a, d, g) = sum(data == d) / length(data);
%             end
% 
% 
%             % OSDS fit
%             if length(select_fov) == 2
%                 data = [OS_DS_direction_total_individual_fit{2, 1, g}{select_fov(1)}, OS_DS_direction_total_individual_fit{2, 1, g}{select_fov(1)} + 180, ...
%                     OS_DS_direction_total_individual_fit{2, 1, g}{select_fov(2)}, OS_DS_direction_total_individual_fit{2, 1, g}{select_fov(2)} + 180];
%             elseif length(select_fov) == 1
%                 data = [OS_DS_direction_total_individual_fit{2, 1, g}{select_fov}, OS_DS_direction_total_individual_fit{2, 1, g}{select_fov} + 180];
%             end
%             tmp_h = polarhistogram(deg2rad(data), 'BinEdges', -1/8*pi : 1/4*pi : 15/8*pi);
%             direction_count_n_animal{8}(a, :, g) = tmp_h.BinCounts / length(data);
% 
% 
%             % OS_only fit
%             if length(select_fov) == 2
%                 data = [OS_only_direction_total_individual_fit{2, 1, g}{select_fov(1)}, OS_only_direction_total_individual_fit{2, 1, g}{select_fov(1)} + 180, ...
%                     OS_only_direction_total_individual_fit{2, 1, g}{select_fov(2)}, OS_only_direction_total_individual_fit{2, 1, g}{select_fov(2)} + 180];
%             elseif length(select_fov) == 1
%                 data = [OS_only_direction_total_individual_fit{2, 1, g}{select_fov}, OS_only_direction_total_individual_fit{2, 1, g}{select_fov} + 180];
%             end
%             tmp_h = polarhistogram(deg2rad(data), 'BinEdges', -1/8*pi : 1/4*pi : 15/8*pi);
%             direction_count_n_animal{9}(a, :, g) = tmp_h.BinCounts / length(data);
% 
% 
%             % OSDS direction fit
%             if length(select_fov) == 2
%                 data = [OS_DS_direction_direction_total_individual_fit{2, 1, g}{select_fov(1)}, OS_DS_direction_direction_total_individual_fit{2, 1, g}{select_fov(2)}];
%             elseif length(select_fov) == 1
%                 data = OS_DS_direction_direction_total_individual_fit{2, 1, g}{select_fov};
%             end
%             tmp_h = polarhistogram(deg2rad(data), 'BinEdges', -1/8*pi : 1/4*pi : 15/8*pi);
%             direction_count_n_animal{10}(a, :, g) = tmp_h.BinCounts / length(data);

            a = a + 1;
        end
    end
end



%% only t2, rose plots or histograms of cell proportions, FOV- or animal-based
mkdir '1_noSmooth\rose_wsem'
g_tag = {'saline', 'gaba'};
t = 2;
for isanimal = [1]
    for cc = 1:4%[3:5, 7:9] 

        for th = 1%:2
            h = figure;
            set(h, 'position', [0 0 1000 1000])

            anova_data_mat = [];
            anova_g_mat{1} = [];
            anova_g_mat{2} = [];
            for g = 1:length(tag)
                switch cc
                    case 1
                        data = DS_direction_total{t, th, g};
                        interval = 1:8;
                        cc_tag = 'DS';
                    case 2
                        proportion_d(g, :) = DS_fit_proportion{th, t}(g, :);
                        cc_tag = 'DS_fit';
                    case 3
                        data = OS_direction_total{t, th, g};
                        data = [data; data + 4];
                        interval = 1:8;
                        cc_tag = 'OS';
                    case 4
                        proportion_d(g, :) = OS_fit_proportion{th, t}(g, :);
                        cc_tag = 'OS_fit';
                    case 5
                        data = OS_DS_direction_total{t, th, g};
                        data = [data; data + 4];
                        interval = 1:8;
                        cc_tag = 'OS_DS';
                    case 6
                        data = OS_DS_direction_direction_total{t, th, g};
                        interval = 1:8;
                        cc_tag = 'OS_DS_direction';
                    case 7
                        data = OS_only_direction_total{t, th, g};
                        data = [data; data + 4];
                        interval = 1:8;
                        cc_tag = 'OS_only';
                    case 8
                        proportion_d(g, :) = fit_OS_DS_proportion{th, t}(g, :);
                        cc_tag = 'OS_DS_fit';
                    case 9
                        proportion_d(g, :) = fit_OS_only_proportion{th, t}(g, :);
                        cc_tag = 'OS_only_fit';
                    case 10
                        proportion_d(g, :) = fit_OS_DS_direction_proportion{th, t}(g, :);
                        cc_tag = 'OS_DS_fit_direction';
                end

                
                if isanimal
                    individual_data{g} = direction_count_n_animal{cc}(:, :, g);
%                     if g == 1
%                         individual_data{g} = individual_data{g}(1:end-1, :);
%                     end
                else
                    individual_data{g} = direction_count_n{cc}(:, :, g);
                end
                notnan_id = find(~isnan(individual_data{g}(:, 1)));
                individual_data{g} = individual_data{g}(notnan_id, :);
                fov_n(cc, g) = length(notnan_id);


                if ismember(cc, [1 3 5:7])
                    [counts_d(g, :), centers_d(g, :)] = hist(data, interval);
                    proportion_d(g, :) = counts_d(g, :) / length(data);
                    if ismember(cc, [3 5 7])
                        proportion_d(g, :) = proportion_d(g, :) * 2;
                    end
                end



                r_vector = ((0 : 1/4 : 2) - 1/2) * pi;
                sem = std(individual_data{g}) ./ sqrt(size(individual_data{g}, 1));
                avg = mean(individual_data{g});
                for d = 1 : 8
                    % use individual mean instead of pooled mean
                    start_p(1, :) = [r_vector(d), avg(d) + sem(d)];
                    start_p(2, :) = [r_vector(d), avg(d) - sem(d)];
                    if d < 8
                        end_p(1, :) = [r_vector(d+1), avg(d+1) + sem(d+1)];
                        end_p(2, :) = [r_vector(d+1), avg(d+1) - sem(d+1)];
                    else
                        end_p(1, :) = [r_vector(1), avg(1) + sem(1)];
                        end_p(2, :) = [r_vector(1), avg(1) - sem(1)];
                    end
                    
%                     start_p(1, :) = [r_vector(d), proportion_d(g, d) + sem(d)];
%                     start_p(2, :) = [r_vector(d), proportion_d(g, d) - sem(d)];
%                     if d < 8
%                         end_p(1, :) = [r_vector(d+1), proportion_d(g, d+1) + sem(d+1)];
%                         end_p(2, :) = [r_vector(d+1), proportion_d(g, d+1) - sem(d+1)];
%                     else
%                         end_p(1, :) = [r_vector(1), proportion_d(g, 1) + sem(1)];
%                         end_p(2, :) = [r_vector(1), proportion_d(g, 1) - sem(1)];
%                     end
                    %                     fill(p(:, 1), p(:, 2), colorVector(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.2); hold on;
                    if ismember(cc, [3:5, 7:9]) && d > 3
                        [x, y] = pol2cart([start_p(1:2, 1); end_p(1:2, 1)], [start_p(1:2, 2); end_p(1:2, 2)]);
                        fill(x([1 2 4 3]), y([1 2 4 3]), colorVector(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.15); hold on;
                    else
                        [x, y] = pol2cart([start_p(1:2, 1); end_p(1:2, 1)], [start_p(1:2, 2); end_p(1:2, 2)]);
                        fill(x([1 2 4 3]), y([1 2 4 3]), colorVector(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on;
                    end
                end


                if g == 1
                    circle_r = max(avg+sem);
%                     circle_r = max(proportion_d(g, :)+sem);
                else
                    circle_r = max(circle_r, max(avg+sem));
%                     circle_r = max(circle_r, max(proportion_d(g, :)+sem));
                end
                avg_p = [(((0 : 1/4 : 2) - 1/2) * pi)', [avg, avg(1)]'];
%                 avg_p = [(((0 : 1/4 : 2) - 1/2) * pi)', [proportion_d(g, :), proportion_d(g, 1)]'];
                [avg_x, avg_y] = pol2cart(avg_p(:, 1), avg_p(:, 2));
                if ismember(cc, [3:5, 7:9])
                    line(avg_x(1:4), avg_y(1:4), 'Color', colorVector(g, :), 'linewidth', 5); hold on
                    line(avg_x(4:end), avg_y(4:end), 'Color', colorVector(g, :), 'linewidth', 2, 'lineStyle', '--'); hold on
                else
                    line(avg_x, avg_y, 'Color', colorVector(g, :), 'linewidth', 5); hold on
                end
                %             polarplot(((0 : 1/4 : 2) - 1/2) * pi, [proportion_d(g, :), proportion_d(g, 1)], 'Color', colorVector(g, :), 'linewidth', 4)
                hold on

                if g == 2
                    myDrawCircle(0, 0, circle_r); hold on
                end
                line([-circle_r, circle_r], [0 0], 'color', [129 129 129]/255, 'lineWidth', 2, 'lineStyle', '--'); hold on
                line([0 0], [-circle_r, circle_r], 'color', [129 129 129]/255, 'lineWidth', 2, 'lineStyle', '--'); hold on
                set(gca,'linew',2, 'FontSize', 14, 'XColor', 'none','YColor','none');
                title([cc_tag, ' th', num2str(th), ' t', num2str(t), ' circle=', num2str(circle_r)])
                axis image


                anova_data_mat = [anova_data_mat; individual_data{g}(:)];
                tmp = repmat(1:8, size(individual_data{g}, 1), 1);
                anova_g_mat{1} = [anova_g_mat{1}; tmp(:)];
                tmp2 = repmat(g_tag(g), size(individual_data{g}, 1)*size(individual_data{g}, 2), 1);
                anova_g_mat{2} = [anova_g_mat{2}; tmp2];

            end
            hold off
            saveas(h, [save_tag, num2str(cc), cc_tag, '_rose_t2only_wsem_animal', num2str(isanimal), '.png'])
            if isanimal
                [d_count_animal_anova_p{cc}{th}, d_count_animal_t{cc}{th}, d_count_animal_stats{cc}{th}] = anovan(anova_data_mat, {anova_g_mat{2}, anova_g_mat{1}}, 'interaction');
                [d_count_animal_c{cc}{th}, m] = multcompare(d_count_animal_stats{cc}{th});
                d_count_animal_anova_p{cc}{th}
            else
                [d_count_anova_p{cc}{th}, d_count_t{cc}{th}, d_count_stats{cc}{th}] = anovan(anova_data_mat, {anova_g_mat{2}, anova_g_mat{1}}, 'interaction');
                [d_count_c{cc}{th}, m] = multcompare(d_count_stats{cc}{th});
                d_count_anova_p{cc}{th}
            end


            for d = 1 : 8
                if isanimal
                    [~, simpleT_animal_pval(cc, th, d)] = ttest2(individual_data{1}(:, d), individual_data{2}(:, d));
                else
                    [~, simpleT_pval(cc, th, d)] = ttest2(individual_data{1}(:, d), individual_data{2}(:, d));
                end
            end



            % cartesian plot for OS
            if ismember(cc, [3:5, 7:9])
                for g = 1:length(tag)
                    switch cc
                        case 3
                            data = OS_direction_total{t, th, g};
                            data = [data; data + 4];
                            interval = 1:8;
                            cc_tag = 'OS';
                        case 4
                            proportion_d(g, :) = OS_fit_proportion{th, t}(g, :);
                            cc_tag = 'OS_fit';
                        case 5
                            data = OS_DS_direction_total{t, th, g};
                            data = [data; data + 4];
                            interval = 1:8;
                            cc_tag = 'OS_DS';
                        case 7
                            data = OS_only_direction_total{t, th, g};
                            data = [data; data + 4];
                            interval = 1:8;
                            cc_tag = 'OS_only';
                        case 8
                            proportion_d(g, :) = fit_OS_DS_proportion{th, t}(g, :);
                            cc_tag = 'OS_DS_fit';
                        case 9
                            proportion_d(g, :) = fit_OS_only_proportion{th, t}(g, :);
                            cc_tag = 'OS_only_fit';
                    end


%                     if isanimal
%                         individual_data{g} = direction_count_n_animal{cc}(:, :, g);
%                     else
%                         individual_data{g} = direction_count_n{cc}(:, :, g);
%                     end
%                     notnan_id = find(~isnan(individual_data{g}(:, 1)));
%                     individual_data{g} = individual_data{g}(notnan_id, :);
%                     fov_n(cc, g) = length(notnan_id);


                    if ismember(cc, [1 3 5:7])
                        [counts_d(g, :), centers_d(g, :)] = hist(data, interval);
                        proportion_d(g, :) = counts_d(g, :) / length(data);
                        if ismember(cc, [3 5 7])
                            proportion_d(g, :) = proportion_d(g, :) * 2;
                        end
                    end


                    r_vector = ((0 : 1/4 : 2) - 1/2) * pi;
                    sem = std(individual_data{g}) ./ sqrt(size(individual_data{g}, 1));
                    avg = mean(individual_data{g});
                    for d = 1 : 3
%                         start_p(1, :) = [d, proportion_d(g, d) + sem(d)];
%                         start_p(2, :) = [d, proportion_d(g, d) - sem(d)];
% 
%                         end_p(1, :) = [d+1, proportion_d(g, d+1) + sem(d+1)];
%                         end_p(2, :) = [d+1, proportion_d(g, d+1) - sem(d+1)];
                        % use individual mean instead of pooled mean
                        fill([d d+1 d+1 d], [avg(d) + sem(d), avg(d+1) + sem(d+1), ...
                            avg(d+1) - sem(d+1), avg(d) - sem(d)], colorVector(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on;
%                         fill([d d+1 d+1 d], [proportion_d(g, d) + sem(d), proportion_d(g, d+1) + sem(d+1), ...
%                             proportion_d(g, d+1) - sem(d+1), proportion_d(g, d) - sem(d)], colorVector(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on;
                    end

                    line(1:4, avg(1:4), 'Color', colorVector(g, :), 'linewidth', 5); hold on
%                     line(1:4, proportion_d(g, 1:4), 'Color', colorVector(g, :), 'linewidth', 5); hold on
                    xlim([0 5])

                    tmp = individual_data{g}(:, 1:4);
                    anova_data_mat = [anova_data_mat; tmp(:)];
                    tmp = repmat(1:4, size(individual_data{g}, 1), 1);
                    anova_g_mat{1} = [anova_g_mat{1}; tmp(:)];
                    tmp2 = repmat(g_tag(g), size(individual_data{g}, 1)*4, 1);
                    anova_g_mat{2} = [anova_g_mat{2}; tmp2];
                end
                hold off
                box off
                set(gca, 'XTick', 0:1:5, 'linew',6, 'FontSize', 22)
                saveas(h, ['1_noSmooth\rose_wsem\', num2str(cc), cc_tag, '_hist_t2only_wsem_animal', num2str(isanimal), '.png'])

                
                if isanimal
                    [o_count_animal_anova_p{cc}{th}, o_count_animal_t{cc}{th}, o_count_animal_stats{cc}{th}] = anovan(anova_data_mat, {anova_g_mat{2}, anova_g_mat{1}}, 'interaction');
                    [o_count_animal_c{cc}{th}, m] = multcompare(o_count_animal_stats{cc}{th});
                    o_count_animal_anova_p{cc}{th}
                else                   
                    [o_count_anova_p{cc}{th}, o_count_t{cc}{th}, o_count_stats{cc}{th}] = anovan(anova_data_mat, {anova_g_mat{2}, anova_g_mat{1}}, 'interaction');
                    [o_count_c{cc}{th}, m] = multcompare(o_count_stats{cc}{th});
                    o_count_anova_p{cc}{th}
                end
            end

        end

    end
end



%% individual median values of DSI/OSI
% col1: DSI, col2: DSI_fit, col3: DSI_global, col4: OSI, col5: OSI_fit, col6: OSI_global
ymax = [0.6 0.6 0.6 0.6 1 0.5];
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
            
            data = squeeze(median_individual{g}(t, cc, :));
            scatter(g * ones(1, length(data)), data, 'MarkerEdgeColor', colorVector(g, :)', 'LineWidth', 2)
            hold on
        end
        
        clear tmp
        tmp{1} = squeeze(median_individual{1}(t, cc, :));
        tmp{2} = squeeze(median_individual{2}(t, cc, :));
        col=@(x)reshape(x,numel(x),1);
        boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
        boxplot2(tmp, 'Width', .6, 'Colors', colors, 'Labels', ftag, 'whisker', 1000)
        box off
        
        title([cc_tag, ' t', num2str(t)])
        xlim([0 length(tag)+1]); ylim([0 1])
        ylabel('proportion')
        set(findobj(gca,'type','line'),'linew',3)
        set(gca,'linew',2, 'FontSize', 14)
        
        [p3(t, cc), hh3(t, cc)] = ttest2(tmp{1}, tmp{2});
        ylim([0 ymax(cc)])
        
    end
    saveas(h, [save_tag, cc_tag, '_median.png'])
    
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
%                 pair_c = pair_counts_OS;
%             case 2
%                 pair_c = pair_counts_fit_OS;
%             case 3
%                 pair_c = pair_counts_DS;
%             case 4
%                 pair_c = pair_counts_fit_DS;
%         end
%         
%         
%         for th = 1 : 2
%             total_pair = length(dist_p{cc, t, th});
%             h = figure; set(h, 'position', [0 0 1000 1000])
%             for dd = 1 : 6
%                 pair_mat{1} = []; pair_mat{2} = [];
%                 for fov_n = 1 : length(pair_c{t, th, 1})
%                     tmp = pair_c{t, th, 1};
%                     if ~isempty(tmp{fov_n}) && sum(tmp{fov_n}{dd}) > 0
%                         pair_mat{1} = [pair_mat{1}; tmp{fov_n}{dd}/sum(tmp{fov_n}{dd})];
%                     end
%                 end
%                 
%                 for fov_n = 1 : length(pair_c{t, th, 2})
%                     tmp = pair_c{t, th, 2};
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
%             saveas(h, [save_tag, 'Spatial_regularity_t', num2str(t), '_', title_tag{cc}, num2str(th), '.png'])
%         end
%     end
% end




%% spatial regularity
% colors = colorVector([1, 2], :);
% title_tag = {'OS', 'OS(f)', 'DS', 'DS(f)'};
% for t = 2
%     for cc = 3 : 4
%         
%         clear pair_mean pair_sem
%         
%         switch cc
%             case 1
%                 pair_c = pair_counts_OS;
%             case 2
%                 pair_c = pair_counts_fit_OS;
%             case 3
%                 pair_c = pair_counts_DS;
%             case 4
%                 pair_c = pair_counts_fit_DS;
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
%                 for fov_n = 1 : length(pair_c{t, th, 1})
%                     tmp = pair_c{t, th, 1};
%                     if ~isempty(tmp{fov_n}) && sum(tmp{fov_n}{dd}) > 0
%                         pair_mat{1} = [pair_mat{1}; tmp{fov_n}{dd}/sum(tmp{fov_n}{dd})];
%                     end
%                 end
%                 
%                 for fov_n = 1 : length(pair_c{t, th, 2})
%                     tmp = pair_c{t, th, 2};
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
%             
%             
%             % plot distance vs difference
%             h = figure;
%             
%             for g = 1 : 2
%                 tmp = pair_c{t, th, g};
%                 for a = 1 : length(tmp)
%                     individual_pair = tmp{a};
%                     for dd = 1 : 5
%                         dist_diff{g}(dd, a) = (individual_pair{dd} * (1 : length(individual_pair{dd}))') / sum(individual_pair{dd});
%                     end
%                 end
%                 
%                 mean_dist_pair(:, g) = mean(dist_diff{g}, 2);
%                 std_dist_pair(:, g) = std(dist_diff{g}, [], 2);
%             end
%             
%             for dd = 1 : 5
%                 [p5{th}(cc, dd), hh5{th}(cc, dd)] = ttest2(dist_diff{1}(dd, :), dist_diff{2}(dd, :));
%             end
%             
%             errorbar(1:5, mean_dist_pair(:, 1), std_dist_pair(:, 1)); hold on
%             errorbar(1:5, mean_dist_pair(:, 2), std_dist_pair(:, 2)); hold on
%             
%         end
%         
%     end
% end
% 





%%
close all
close all hidden
fn_tmp = {'noSmoothing', 'smoothing'};
save([save_tag, 'stats_', fn_tmp{is_smooth+1}, '_th', num2str(ii), '_FOV.mat'], 'hh1', 'hh2', 'hh3', 'o_count_animal_anova_p', 'o_count_animal_t', 'o_count_animal_stats', 'o_count_animal_c', ...
    'd_count_animal_anova_p', 'd_count_animal_t', 'd_count_animal_stats', 'd_count_animal_c', 'simpleT_animal_pval')

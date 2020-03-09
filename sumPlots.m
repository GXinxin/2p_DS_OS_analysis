%% Plot dist vs similarity for all cell pairs
%     colorAdd = [228 180 34]/255;
colorAdd = [1 0 0];
titleText = {'saline', 'gabazine'};
%     cc_tag = {'OS', 'OS(f)', 'DS', 'DS(f)'};
%     max_value = [2 90 90 4 180];
cc_tag = {'OS', 'OS(waxis)', 'OS(wdir)', 'OSonly', 'OSonly(waxis)', 'OS_only(wdir)', 'OSDS', 'OSDS(waxis)', 'OS_DS(wdir)', ...
    'OSfit', 'OSonlyFit', 'OSDSFit', 'DS', 'weightedDs', 'DSfit', 'gOSI', 'gDSI', 'all'};
max_value = [2 90 90 2 90 90 2 90 90 90 90 90 4 180 180 2 4 4];
max_bin = 8;
clear tmp

b_real = cell(4, 2, 2);
bint_real = cell(4, 2, 2);
r_real = cell(4, 2, 2);
rint_real = cell(4, 2, 2);
stats = cell(4, 2, 2);
b_test = cell(4, 2);
bint_test = cell(4, 2);
r_test = cell(4, 2);
rint_test = cell(4, 2);
stats_test = cell(4, 2);

for case_i = [1 4 7 16]
    bin_centers{case_i} = 0:1:2;
end
for case_i = [13 17 18]
    bin_centers{case_i} = 0:1:4;
end
for case_i = [2 3 5 6 8 9 10:12]
    bin_centers{case_i} = 5 : 10 : 85;
end
for case_i = 14:15
    bin_centers{case_i} = 10 : 20 : 180;
end

clear case_i

t = 2;
for case_i = 16:18%[1 3 4 6 7 9:17]
    for th = 1 %: 2
        h = figure; set(h, 'position', [0 0 1200 500])
        for g = 1 : 2
            subplot(1, 2, g)
            scatter(total_pdist{th, g}{case_i}(:, 1) * pixel_res, total_pdist{th, g}{case_i}(:, 2), '.',...
                'MarkerFaceColor', colorVector(g, :), 'MarkerEdgeColor', colorVector(g, :)); hold on
            title([titleText{g}, ' ' cc_tag{case_i}, ' th', num2str(th), ' t', num2str(t), ' ', num2str(size(total_pdist{th, g}{case_i}, 1))])
            
            line([0 600*pixel_res], [max_value(case_i)/2, max_value(case_i)/2], 'color', 'k', 'LineStyle','--'); hold on
            
            dist_bin = ceil(total_pdist{th, g}{case_i}(:, 1) / bin_pixel);
            for b = 1 : max_bin % max(dist_bin)
                p_id = find(dist_bin == b);
                avg_diff(case_i, b) = mean(total_pdist{th, g}{case_i}(p_id, 2));
                sem_diff(case_i, b) = std(total_pdist{th, g}{case_i}(p_id, 2))/sqrt(length(total_pdist{th, g}{case_i}(p_id, 2)));
                x = [(b-1)*bin_um, b*bin_um];
                y = [avg_diff(case_i, b), avg_diff(case_i, b)];
                curve1 = y + sem_diff(case_i, b);
                curve2 = y - sem_diff(case_i, b);
                x2 = [x, fliplr(x)];
                inBetween = [curve1, fliplr(curve2)];
                fill(x2, inBetween, colorAdd, 'FaceAlpha', 0.2, 'EdgeColor', colorAdd); hold on;
                line(x, y, 'color', colorAdd, 'lineWidth', 1.5); hold on
                if b > 1
                    line([(b-1)*bin_um-1, (b-1)*bin_um+1], [avg_diff(case_i, b-1), avg_diff(case_i, b)], 'color', colorAdd, 'lineWidth', 1); hold on
                end
                
                % individual points within each bin
                bin_value{case_i}{g}{b} = total_pdist{th, g}{case_i}(p_id, 2);
                hist_bin_value{case_i}{g}(b, :) = hist(bin_value{case_i}{g}{b}, bin_centers{case_i});
            end
            
            X = [total_pdist{th, g}{case_i}(:, 1), ones(length(total_pdist{th, g}{case_i}(:, 1)), 1)];
            y = total_pdist{th, g}{case_i}(:, 2);
            [b_real{case_i, th, g}, bint_real{case_i, th, g}, r_real{case_i, th, g}, rint_real{case_i, th, g}, stats{case_i, th, g}] = regress(y, X);
            x = 0 : (500 * pixel_res);
            y = b_real{case_i, th, g}(1) * x + b_real{case_i, th, g}(2);
            plot(x, y, 'color', 'k', 'LineWidth', 3);
            set(gca,'linew',2, 'FontSize', 10)
            xlim([0 500]);
            
            if g == 1
                y = [total_pdist{th, 1}{case_i}(:, 2);  total_pdist{th, 2}{case_i}(:, 2)];
                X1 = [total_pdist{th, 1}{case_i}(:, 1); total_pdist{th, 2}{case_i}(:, 1)];
                X2 = [ones(length(total_pdist{th, 1}{case_i}(:, 1)), 1); 2*ones(length(total_pdist{th, 2}{case_i}(:, 1)), 1)];
                X = [X1, X2, X1.*X2, ones(length(X1), 1)];
                [b_test{case_i, th}, bint_test{case_i, th}, r_test{case_i, th}, rint_test{case_i, th}, stats_test{case_i, th}] = regress(y, X);
            end
        end
        saveas(h, [save_tag, 'DistVsDiff', titleText{g}, '_', cc_tag{case_i}, '_th', num2str(th), '_t', num2str(t), '.png'])
        
        
        % histogram of pref_diff in several distance bins
        facevalue = [0.3 0];
        h = figure; set(h, 'position', [500 700 1600 500])
        for b = 1 : max_bin
            subplot(2, 4, b)
            for g = 1 : 2
                bar(bin_centers{case_i}, hist_bin_value{case_i}{g}(b, :)/sum(hist_bin_value{case_i}{g}(b, :)), 'barWidth', 0.9, 'FaceColor', colorVector(g, :),...
                    'FaceAlpha', facevalue(g), 'EdgeColor', colorVector_line(g, :), 'LineWidth', 1.5); hold on
            end
            set(gca,'FontSize', 7.5);
            title([cc_tag{case_i}, ' th', num2str(th)])
        end
        saveas(h, [save_tag, 'DistVsDiff_corrected', titleText{g}, '_', cc_tag{case_i}, '_th', num2str(th), '_hist.png'])
        
        
        % histogram of pref_diff in several distance bins
        h = figure; set(h, 'position', [150 700 1500 500])
        for b = 1 : max_bin
            subplot(2, 4, b)
            for g = 1 : 2
                if ~isempty(bin_value{case_i}{g}{b})
                    ch = cdfplot(bin_value{case_i}{g}{b});
                    ch.Color = colorVector(g, :);
                    ch.LineWidth = 3; hold on
                    if g == 2
                        [~, ks_p{case_i}(th, b)] = kstest2(bin_value{case_i}{1}{b}, bin_value{case_i}{2}{b});
                    end
                end
            end
            title([cc_tag{case_i}, ' p=', num2str(ks_p{case_i}(th, b)), ' th', num2str(th)])
        end
        saveas(h, [save_tag, 'DistVsDiff', titleText{g}, '_', cc_tag{case_i}, '_th', num2str(th), '_cdf.png'])
    end
end

save([save_tag, 'fittedSlope.mat'], 'b_real', 'bint_real', 'r_real', 'rint_real', 'stats', 'b_test', 'bint_test', 'r_test', 'rint_test', 'stats_test')





%% Plot dist vs similarity for all cell pairs: correlation coeffient
colorAdd = [1 0 0];
titleText = {'saline', 'gabazine'};
cc_tag = {'OS(axis)', 'OS(dir)', 'OSonly(axis)', 'OSonly(dir)', 'OSDS(axis)', 'OSDS(dir)', 'DS', 'all', 'gOSI', 'gDSI'};
max_value = ones(1, length(cc_tag));
max_bin = 8;
clear tmp
load('E:\Lab\Data\2p\tra2b\summary\1_noSmooth\p15_tra2b_fittedSlope_corr.mat')

b_real = cell(4, 2, 2);
bint_real = cell(4, 2, 2);
r_real = cell(4, 2, 2);
rint_real = cell(4, 2, 2);
stats = cell(4, 2, 2);
b_test = cell(4, 2);
bint_test = cell(4, 2);
r_test = cell(4, 2);
rint_test = cell(4, 2);
stats_test = cell(4, 2);


for case_i = 1:10
    bin_centers_corr{case_i} = -0.9:0.2:1;
end


t = 2;
for case_i = 8: 10
    for th = 1 %: 2
        h = figure; set(h, 'position', [0 0 1200 500])
        for g = 1 : 2
            subplot(1, 2, g)
            scatter(total_corr{th, g}{case_i}(:, 1) * pixel_res, total_corr{th, g}{case_i}(:, 2), '.',...
                'MarkerFaceColor', colorVector(g, :), 'MarkerEdgeColor', colorVector(g, :)); hold on
            title([titleText{g}, ' ' cc_tag{case_i}, ' th', num2str(th), ' t', num2str(t), ' ', num2str(size(total_corr{th, g}{case_i}, 1))])
            
            line([0 600*pixel_res], [0,0], 'color', 'k', 'LineStyle','--'); hold on
            
            dist_bin = ceil(total_corr{th, g}{case_i}(:, 1) / bin_pixel);
            for b = 1 : max_bin % max(dist_bin)
                p_id = find(dist_bin == b);
                avg_diff(case_i, b) = mean(total_corr{th, g}{case_i}(p_id, 2));
                sem_diff(case_i, b) = std(total_corr{th, g}{case_i}(p_id, 2))/sqrt(length(total_corr{th, g}{case_i}(p_id, 2)));
                x = [(b-1)*bin_um, b*bin_um];
                y = [avg_diff(case_i, b), avg_diff(case_i, b)];
                curve1 = y + sem_diff(case_i, b);
                curve2 = y - sem_diff(case_i, b);
                x2 = [x, fliplr(x)];
                inBetween = [curve1, fliplr(curve2)];
                fill(x2, inBetween, colorAdd, 'FaceAlpha', 0.2, 'EdgeColor', colorAdd); hold on;
                line(x, y, 'color', colorAdd, 'lineWidth', 1.5); hold on
                if b > 1
                    line([(b-1)*bin_um-1, (b-1)*bin_um+1], [avg_diff(case_i, b-1), avg_diff(case_i, b)], 'color', colorAdd, 'lineWidth', 1.5); hold on
                end
                
                % individual points within each bin
                bin_value_corr{case_i}{g}{b} = total_corr{th, g}{case_i}(p_id, 2);
                hist_bin_value_corr{case_i}{g}(b, :) = hist(bin_value_corr{case_i}{g}{b}, bin_centers_corr{case_i});
            end
            
            X = [total_corr{th, g}{case_i}(:, 1), ones(length(total_corr{th, g}{case_i}(:, 1)), 1)];
            y = total_corr{th, g}{case_i}(:, 2);
            [b_real{case_i, th, g}, bint_real{case_i, th, g}, r_real{case_i, th, g}, rint_real{case_i, th, g}, stats{case_i, th, g}] = regress(y, X);
            x = 0 : (600 * pixel_res);
            y = b_real{case_i, th, g}(1) * x + b_real{case_i, th, g}(2);
            plot(x, y, 'color', 'k', 'LineWidth', 2);
            xlim([0 500]);
            ylim([-1 1]);
            set(gca,'linew',2, 'FontSize', 10)
            
            if g == 1
                y = [total_corr{th, 1}{case_i}(:, 2);  total_corr{th, 2}{case_i}(:, 2)];
                X1 = [total_corr{th, 1}{case_i}(:, 1); total_corr{th, 2}{case_i}(:, 1)];
                X2 = [ones(length(total_corr{th, 1}{case_i}(:, 1)), 1); 2*ones(length(total_corr{th, 2}{case_i}(:, 1)), 1)];
                X = [X1, X2, X1.*X2, ones(length(X1), 1)];
                [b_test{case_i, th}, bint_test{case_i, th}, r_test{case_i, th}, rint_test{case_i, th}, stats_test{case_i, th}] = regress(y, X);
            end
        end
        saveas(h, [save_tag, 'DistVsCorr_correct', titleText{g}, '_', cc_tag{case_i}, '_th', num2str(th), '_t', num2str(t), '.png'])
        
        
        % histogram of pref_diff in several distance bins
        facevalue = [0.3 0];
        h = figure; set(h, 'position', [500 700 1600 500])
        for b = 1 : max_bin
            subplot(2, 4, b)
            for g = 1 : 2
                bar(bin_centers_corr{case_i}, hist_bin_value_corr{case_i}{g}(b, :)/sum(hist_bin_value_corr{case_i}{g}(b, :)), 'barWidth', 0.9, 'FaceColor', colorVector(g, :),...
                    'FaceAlpha', facevalue(g), 'EdgeColor', colorVector_line(g, :), 'LineWidth', 1.5); hold on
            end
            set(gca,'FontSize', 7.5);
            title([cc_tag{case_i}, ' th', num2str(th)])
        end
        saveas(h, [save_tag, 'DistVsCorr', titleText{g}, '_', cc_tag{case_i}, '_th', num2str(th), '_hist.png'])
        
        
        % cdf of pref_diff
        h = figure; set(h, 'position', [150 700 1500 500])
        for b = 1 : max_bin
            subplot(2, 4, b)
            for g = 1 : 2
                if ~isempty(bin_value_corr{case_i}{g}{b})
                    ch = cdfplot(bin_value_corr{case_i}{g}{b});
                    ch.Color = colorVector(g, :);
                    ch.LineWidth = 3; hold on
                    if g == 2
                        [~, ks_p_corr{case_i}(th, b)] = kstest2(bin_value_corr{case_i}{1}{b}, bin_value_corr{case_i}{2}{b});
                    end
                end
            end
            title([cc_tag{case_i}, ' p=', num2str(ks_p_corr{case_i}(th, b)), ' th', num2str(th)])
        end
        saveas(h, [save_tag, 'DistVsCorr', titleText{g}, '_', cc_tag{case_i}, '_th', num2str(th), '_cdf.png'])
    end
end

save([save_tag, 'fittedSlope_corr.mat'], 'b_real', 'bint_real', 'r_real', 'rint_real', 'stats', 'b_test', 'bint_test', 'r_test', 'rint_test', 'stats_test')




%% concatenate data from animals for circular variance and angle proportion
for g = 1 : 2
    animal_c = 1;
    for ff = 1 : size(is_include{g}, 1)
        if sum(is_include{g}(ff, :)) > 0
            for case_i = 1 : 10
                if sum(is_include{g}(ff, :)) == 2
                    if ismember(case_i, [1 3 5 6 7])
                        tmp_concat{g}{animal_c} = [tmp_animal{g}{case_i}{ff, 1}; tmp_animal{g}{case_i}{ff, 2}];
                        tmp_concat2{g}{animal_c} = [tmp_animal2{g}{case_i}{ff, 1}; tmp_animal2{g}{case_i}{ff, 2}];
                    else
                        tmp_concat{g}{animal_c} = [tmp_animal{g}{case_i}{ff, 1}; tmp_animal{g}{case_i}{ff, 2}];
                        tmp_concat2{g}{animal_c} = [tmp_animal2{g}{case_i}{ff, 1}, tmp_animal2{g}{case_i}{ff, 2}];
                    end
                elseif sum(is_include{g}(ff, :)) == 1
                    id = find(is_include{g}(ff, :) == 1);
                    tmp_concat{g}{animal_c} = tmp_animal{g}{case_i}{ff, id};
                    tmp_concat2{g}{animal_c} = tmp_animal2{g}{case_i}{ff, id};
                end
                
                if ~isempty(tmp_concat{g}{animal_c})
                    fov_circ_var_animal{g}(animal_c, case_i) = circ_var(tmp_concat{g}{animal_c});
                end
                
                
                if ~isempty(tmp_concat{g}{animal_c})
                    if ismember(case_i, [1 5 6 11])
                        within_n = sum(tmp_concat2{g}{animal_c} == 1);
                    elseif ismember(case_i, [4 10])
                        within_n = sum((tmp_concat2{g}{animal_c} < 202.5) & (tmp_concat2{g}{animal_c} > 157.5));
                    elseif ismember(case_i, [2 8 9])
                        within_n = sum((tmp_concat2{g}{animal_c} < 22.5) + (tmp_concat2{g}{animal_c} > 157.5));
                    else
                        within_n = sum(tmp_concat2{g}{animal_c} == 5);
                    end
                    angle_proportion_animal{g}(animal_c, case_i) = within_n / length(tmp_concat2{g}{animal_c});
                end
            end
            animal_c = animal_c + 1;
        end
    end
end






%% Plot circular variance difference between groups
titleText = {'saline', 'gabazine'};
cc_tag = {'OS', 'fitOS', 'DS', 'fitDS', 'OSDS', 'OSonly', 'OSDSdirection', 'fitOSDS', 'fitOSonly', 'fitOSDSdirection'};
tag = {'FOV', 'animal'};
t = 2;
clear tmp
for is_animal = [0 1]
    
    for case_i = 1 : 10
        for th = 1 %: 2
            h = figure;
            for g = 1 : 2
                if is_animal
%                     if g == 2
%                         tmp{g} = fov_circ_var_animal{g}([1:4, 6:8], case_i);
%                     else
                        tmp{g} = fov_circ_var_animal{g}(:, case_i);
%                     end
                else
%                     if g == 2
%                         tmp{g} = fov_circ_var{t, th, g}([1:7, 9:14], case_i);
%                     else
                        tmp{g} = fov_circ_var{t, th, g}(:, case_i);
%                     end
                end
                %                 tmp{g} = tmp{g}(tmp{g} > 0);
                scatter(g * ones(1, length(tmp{g})), tmp{g}, 'MarkerEdgeColor', colorVector(g, :)', 'LineWidth', 2)
                hold on
            end
            
            col=@(x)reshape(x,numel(x),1);
            boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
            boxplot2(tmp, 'Width', .6, 'Colors', colorVector, 'Labels', titleText, 'whisker', 1000)
            
            [p00{is_animal+1}(case_i, th), hh00{is_animal+1}(case_i, th)] = ttest2(tmp{1}, tmp{2});
            
            title([cc_tag{case_i}, ' th', num2str(th), ' ', num2str(hh00{is_animal+1}(case_i, th)), tag{is_animal+1}])
            ylabel('circ var')
            ylim([0 1])
            set(findobj(gca,'type','line'),'linew',3)
            set(gca,'linew',4, 'FontSize', 14)
            box off
            
            saveas(h, [save_tag, cc_tag{case_i}, 'CircVar_t', num2str(t), '_th', num2str(th),'_animal', num2str(is_animal), '.png'])
            
        end
    end
end


%% Plot angle proportion in the opposite to dominant direction
titleText = {'saline', 'gabazine'};
cc_tag = {'OS', 'fitOS', 'DS', 'fitDS', 'OSDS', 'OSonly', 'OSDSdirection', 'fitOSDS', 'fitOSonly', 'fitOSDSdirection'};
tag = {'FOV', 'animal'};
t = 2;
clear tmp
for is_animal = [0 1]
    for case_i = 1 : 10
        for th = 1 %: 2
            h = figure;
            for g = 1 : 2
                if is_animal
%                     if g == 2
%                         tmp{g} = angle_proportion_animal{g}([1:4, 6:8], case_i);
%                     else
                        tmp{g} = angle_proportion_animal{g}(:, case_i);
%                     end
                else
%                     if g == 2
%                         tmp{g} = angle_proportion{t, th, g}([1:7, 9:14], case_i);
%                     else
                        tmp{g} = angle_proportion{t, th, g}(:, case_i);
%                     end
                end
                %                 tmp{g} = tmp{g}(find(tmp{g} > 0));
                scatter(g * ones(1, length(tmp{g})), tmp{g}, 'MarkerEdgeColor', colorVector(g, :)', 'LineWidth', 2)
                hold on
            end
            
            col=@(x)reshape(x,numel(x),1);
            boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
            boxplot2(tmp, 'Width', .6, 'Colors', colorVector, 'Labels', titleText, 'whisker', 1000)
            
            [p01{is_animal+1}(case_i, th), hh01{is_animal+1}(case_i, th)] = ttest2(tmp{1}, tmp{2});
            
            title([cc_tag{case_i}, ' th', num2str(th), ' ', num2str(hh01{is_animal+1}(case_i, th)), tag{is_animal+1}])
            ylabel('angleProportion')
            ylim([0 1])
            set(findobj(gca,'type','line'),'linew',3)
            set(gca,'linew',4, 'FontSize', 14)
            box off
            
            saveas(h, [save_tag, cc_tag{case_i}, 'angleProportion_t', num2str(t), '_th', num2str(th),'_animal', num2str(is_animal), '.png'])
            
        end
    end
end



%% DSI, OSI, time2peak/time2slope
% titleText = {'glutamate blocker', 'cholinergic blocker', 'saline'};
% titleText = {'glutamate blocker', 'saline'};
titleText = {'saline', 'gabazine'};

clear counts centers proportion
for cc = 1 : 8
    
    h = figure; set(h, 'position', [0 0 1000, 1000])
    count = 1;
    for t = 1 : 3
        for g = [1 2]
            
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
        [cc mean(pre_data) std(pre_data)/sqrt(length(pre_data)) mean(data) std(data)/sqrt(length(data))]
    end
    saveas(h, [save_tag, [cc_tag, '_hist.png']])
    
end



%% combined histograms

titleText = {'saline', 'gabazine'};
range = [0.2, 0.2, 0.4, 0.3, 0.3, 0.5];

clear counts centers proportion

for cc = [3 6]%1 : 6
    
    h = figure; set(h, 'position', [0 0 800, 800])
    count = 1;
    t = 2;
    
    for g = [1 2]
        
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
        end
        
        
        median_data = median_total{g}(t, cc);
        if g == 1
            pre_data = data;
            median_data1 = median_data;
        end
        
        %             subplot(3, 2, count)
        if cc > 6
            [counts(g, :), centers(g, :)] = hist(data, 20);
        else
            %                 [counts(g, :), centers(g, :)] = hist(data, 0.025:0.05:1);
            [counts(g, :), centers(g, :)] = hist(data, 0.025:0.05:0.8);
        end
        proportion(g, :) = counts(g, :) / length(data);
        if g == 1
            facevalue = 0.3;
        else
            facevalue = 0.2;
        end
        bar(centers(g, :), proportion(g, :), 'barWidth', 0.9, 'FaceColor', colorVector(g, :),...
            'FaceAlpha', facevalue, 'EdgeColor', colorVector_line(g, :), 'LineWidth', 1.5)
        line(median_data * [1 1], [0 range(cc)], 'lineStyle', '- -', 'color', colorVector(g, :), 'lineWidth', 3)
        if cc > 6
            xlim([0 150])
        else
            xlim([0 0.8])
        end
        
        ylim([0 0.3])
        ylabel('proportion', 'FontSize', 12, 'FontWeight', 'Bold')
        xlabel(cc_tag, 'FontSize', 12, 'FontWeight', 'Bold')
        title([' median=', num2str(median_data1), '  ', num2str(median_data)], 'FontSize', 15)
        set(gca,'linew',3, 'FontSize', 22)
        count = count + 1;
        box off
        hold on
    end
    
    [p1(t, cc), hh1(t, cc)] = ttest2(pre_data, data);
    saveas(h, [save_tag, [cc_tag, '_combined_hist.png']])
    
end




%% combined cdf

titleText = {'saline', 'gabazine'};
range = [0.2, 0.2, 0.4, 0.3, 0.3, 0.5];

clear counts centers proportion

for cc = 1 : 6
    
    h = figure; set(h, 'position', [0 0 500, 500])
    count = 1;
    t = 2;
    
    for g = [1 2]
        
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
        end
        
        
        median_data = median_total{g}(t, cc);
        if g == 1
            pre_data = data;
        end
        
        ch = cdfplot(data);
        ch.Color = colorVector(g, :);
        ch.LineWidth = 3;
        %             bar(centers(g, :), proportion(g, :), 'barWidth', 0.9, 'FaceColor', colorVector(g, :),...
        %                 'FaceAlpha', facevalue, 'EdgeColor', colorVector_line(g, :), 'LineWidth', 1.5)
        %             line(median_data * [1 1], [0 range(cc)], 'lineStyle', '- -', 'color', colorVector(g, :), 'lineWidth', 3)
        if cc > 6
            xlim([0 150])
        else
            xlim([0 1])
        end
        hold on
    end
    
    ylabel('proportion', 'FontSize', 12, 'FontWeight', 'Bold')
    xlabel(cc_tag, 'FontSize', 12, 'FontWeight', 'Bold')
    [p0(cc), hh0(cc)] = kstest2(pre_data, data);
    title(num2str(hh0(cc)), 'FontSize', 15)
    set(gca,'linew',3, 'FontSize', 22)
    box off; grid off
    
    
    saveas(h, [save_tag, [cc_tag, '_combined_cdf.png']])
    
end

close all


perDirectionAnalysis;



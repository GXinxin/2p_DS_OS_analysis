% plots for DSI and OSI quantifications of each FOV



% get spatial location (center point) in the FOV for each cell
tmp = dir(fullfile('ROIs_DS*'));
if ~isempty(tmp)
    cd(tmp(1).name)
    fn = dir(fullfile('*.mldatx'));
    load([fn(1).name, '\all_vars.mat']);
else
    tmp = dir(fullfile('*proc*.mat'));
    load(tmp(1).name)
    for i = 1 : length(dat.stat)
        iscell(i) = dat.stat(i).iscell;
        ROIy(i) = dat.stat(i).med(1);
        ROIx(i) = dat.stat(i).med(2);
    end
    ROIy = ROIy(iscell > 0); ROIy = ROIy(active_id);
    ROIx = ROIx(iscell > 0); ROIx = ROIx(active_id);
end



% plot skewness distribution
h = figure; hist(skew_mat, 0.5:1:15); xlabel('skewness'); ylabel('cell no')
line(median(skew_mat) * [1 1], [0 50], 'color', 'r', 'lineWidth', 2)
title(['median = ', num2str(median(skew_mat))]); xlim([0 16])
saveas(h, [save_folder, save_tag, 'skewness.png'])



% plot DSI and OSI related quantification for each FOV
for t = 1 : 3
    
    colorVector = copper(4);
    
    
    % plot OSI histograms
    h = figure; set(h, 'position', [200 0 500 900])
    subplot(3, 1, 1)
    [counts, centers] = hist(OSI(id1{t}, t), 0.05:0.05:1);
    proportion = counts / length(OSI(id1{t}, t));
    bar(centers, proportion, 'barWidth', 0.9)
    line(median_values(t, 4) * [1 1], [0 max(proportion)], 'color', 'r', 'lineWidth', 2)
    title(['median = ', num2str(median_values(t, 4))])
    xlim([0 1])
    ylim([0 .5])
    ylabel('proportion', 'FontSize', 12, 'FontWeight', 'Bold')
    xlabel('OSI', 'FontSize', 12, 'FontWeight', 'Bold')
    box off
    
    subplot(3, 1, 2)
    [counts, centers] = hist(OSI_fit(id2{t}, t), 0.05:0.05:1);
    proportion = counts / length(OSI_fit(id2{t}, t));
    bar(centers, proportion, 'barWidth', 0.9)
    line(median_values(t, 5) * [1 1], [0 max(proportion)], 'color', 'r', 'lineWidth', 2)
    title(['median = ', num2str(median_values(t, 5))])
    xlim([0 1])
    ylim([0 .5])
    ylabel('proportion', 'FontSize', 12, 'FontWeight', 'Bold')
    xlabel('OSI', 'FontSize', 12, 'FontWeight', 'Bold')
    box off
    
    subplot(3, 1, 3)
    [counts, centers] = hist(OSI_global(id1{t}, t), 0.05:0.05:1);
    proportion = counts / length(OSI_global(id1{t}, t));
    bar(centers, proportion, 'barWidth', 0.9)
    line(median_values(t, 6) * [1 1], [0 max(proportion)], 'color', 'r', 'lineWidth', 2)
    title(['median = ', num2str(median_values(t, 6))])
    xlim([0 1])
    ylim([0 .5])
    ylabel('proportion', 'FontSize', 12, 'FontWeight', 'Bold')
    xlabel('global OSI', 'FontSize', 12, 'FontWeight', 'Bold')
    box off
    saveas(h, [save_folder, save_tag, 'OSI_hist_p', num2str(t), 'fit0.png'])
    
    
    
    % plot preferred orientation
    for th = 1 : 2
        h = figure;
        subplot(2, 1, 1)
        [counts_d, centers_d] = hist((OS_direction{t}{th}-1)*45, 10:10:180);
        proportion_d = counts_d / length(OS_direction{t}{th});
        bar(centers_d, proportion_d, 'barWidth', 1)
        box off
        set(gca,'linew',2, 'FontSize', 15)
        xlim([0 180]); title(['th', num2str(th), ' pref orientation'])
        
        subplot(2, 1, 2)
        [counts_d, centers_d] = hist(fit_OS_direction{t}{th}, 10:10:180);
        proportion_d = counts_d / length(fit_OS_direction{t}{th});
        bar(centers_d, proportion_d, 'barWidth', 1)
        box off
        set(gca,'linew',2, 'FontSize', 15)
        xlim([0 180]); title(['th', num2str(th), ' pref orientation (fitted)'])
        saveas(h, [save_folder, save_tag, 'distribution_orientation_p', num2str(t), '_th', num2str(th), '.png'])
        
        
        % plot preferred orientation as rose
        h = figure; set(h, 'position', [300 100 1000 500])
        subplot(1, 2, 1)
        [counts_d, centers_d] = hist(OS_direction{t}{th}, 1:8);
        proportion_d = counts_d / length(OS_direction{t}{th});
        polarplot(((0 : 1/4 : 2) - 1/2) * pi, [proportion_d, proportion_d(1)], 'linewidth', 4)
        thetalim([-90 90])
        title('OS')
        set(gca,'linew',2, 'FontSize', 15)

        subplot(1, 2, 2)
        [counts_d, centers_d] = hist(fit_OS_direction{t}{th}, 0:45:315);
        proportion_d = counts_d / length(fit_OS_direction{t}{th});
        polarplot(((0 : 1/4 : 2) - 1/2) * pi, [proportion_d, proportion_d(1)], 'linewidth', 4)
        thetalim([-90 90])
        title('OS bimodal fitting')
        set(gca,'linew',2, 'FontSize', 15)
        saveas(h, [save_folder, save_tag, 'OS_rose_p', num2str(t), '_th', num2str(th), '.png'])
        
        
        
        % plot spatial location and preferred orientation together
        colorRange = hsv(181);
        h = figure; set(h, 'position', [0 0 1000 350])
        subplot(1, 2, 1)
        select_id = id1{t}(OS_id{t}{th});
        for c = 1 : length(OS_direction{t}{th})
            scatter(ROIx(select_id(c)), ROIy(select_id(c)), 'MarkerEdgeColor', colorRange(1+uint16((OS_direction{t}{th}(c)-1)*45), :), ...
                'MarkerFaceColor', colorRange(1+uint16((OS_direction{t}{th}(c)-1)*45), :)); hold on
        end
        xlim([0 512]); ylim([0 512]); %axis image
        colormap hsv; colorbar('Ticks',0:0.25:1, 'TickLabels',{'0','45','90','135','180'}); 
        title(['th', num2str(th), ' pref orientation'])
        
        subplot(1, 2, 2)
        select_id = id2{t}(fit_OS_id{t}{th});
        for c = 1 : length(fit_OS_direction{t}{th})
            scatter(ROIx(select_id(c)), ROIy(select_id(c)), 'MarkerEdgeColor', colorRange(1+uint16(fit_OS_direction{t}{th}(c)), :), ...
                'MarkerFaceColor', colorRange(1+uint16(fit_OS_direction{t}{th}(c)), :)); hold on
        end
        xlim([0 512]); ylim([0 512]); %axis image
        colormap hsv; colorbar('Ticks',0:0.25:1, 'TickLabels',{'0','45','90','135','180'}); 
        title(['th', num2str(th), ' pref orientation (fitted)'])
        saveas(h, [save_folder, save_tag, 'spatial_orientation_p', num2str(t), '_th', num2str(th), '.png'])
    end
    
    
    
    if t > 1
        % plot DSI histograms
        h = figure; set(h, 'position', [200 0 500 900])
        subplot(3, 1, 1)
        [counts, centers] = hist(DSI(id1{t}, t), 0.05:0.05:1);
        proportion = counts / length(DSI(id1{t}, t));
        bar(centers, proportion, 'barWidth', 0.9)
        line(median_values(t, 1) * [1 1], [0 max(proportion)], 'color', 'r', 'lineWidth', 2)
        title(['median = ', num2str(median_values(t, 1))])
        xlim([0 1])
        ylim([0 0.6])
        ylabel('proportion', 'FontSize', 12, 'FontWeight', 'Bold')
        xlabel('DSI', 'FontSize', 12, 'FontWeight', 'Bold')
        box off
        
        subplot(3, 1, 2)
        [counts, centers] = hist(DSI_fit(id2{t}, t), 0.05:0.05:1);
        proportion = counts / length(DSI_fit(id2{t}, t));
        bar(centers, proportion, 'barWidth', 0.9)
        line(median_values(t, 2) * [1 1], [0 max(proportion)], 'color', 'r', 'lineWidth', 2)
        title(['median = ', num2str(median_values(t, 2))])
        xlim([0 1])
        ylim([0 0.6])
        ylabel('proportion', 'FontSize', 12, 'FontWeight', 'Bold')
        xlabel('DSI', 'FontSize', 12, 'FontWeight', 'Bold')
        box off
        
        subplot(3, 1, 3)
        [counts, centers] = hist(DSI_global(id1{t}, t), 0.05:0.05:1);
        proportion = counts / length(DSI_global(id1{t}, t));
        bar(centers, proportion, 'barWidth', 0.9)
        line(median_values(t, 3) * [1 1], [0 max(proportion)], 'color', 'r', 'lineWidth', 2)
        title(['median = ', num2str(median_values(t, 3))])
        xlim([0 1])
        ylim([0 0.6])
        ylabel('proportion', 'FontSize', 12, 'FontWeight', 'Bold')
        xlabel('global DSI', 'FontSize', 12, 'FontWeight', 'Bold')
        box off
        saveas(h, [save_folder, save_tag, 'DSI_hist_p', num2str(t), 'fit0.png'])
        
        
        
        for th = 1 : 2
            % plot preferred direction
            h = figure; set(h, 'position', [300 100 1000 500])
            subplot(1, 2, 1)
            [counts_d, centers_d] = hist(DS_direction{t}{th}, 1:8);
            proportion_d = counts_d / length(DS_direction{t}{th});
            polarplot(((0 : 1/4 : 2) - 1/2) * pi, [proportion_d, proportion_d(1)], 'linewidth', 4)
            title('no fitting')
            set(gca,'linew',2, 'FontSize', 15)
            
            subplot(1, 2, 2)
            [counts_d, centers_d] = hist(fit_DS_direction{t}{th}, 0:45:315);
            proportion_d = counts_d / length(fit_DS_direction{t}{th});
            polarplot(((0 : 1/4 : 2) - 1/2) * pi, [proportion_d, proportion_d(1)], 'linewidth', 4)
            title('bimodal fitting')
            set(gca,'linew',2, 'FontSize', 15)
            saveas(h, [save_folder, save_tag, 'distribution_p', num2str(t), '_th', num2str(th), '.png'])
            
            
            % plot spatial location and preferred direction together
            colorRange = hsv(8);
            h = figure; set(h, 'position', [0 0 1000 350])
            subplot(1, 2, 1)
            select_id = id1{t}(DS_id{t}{th});
            for c = 1 : length(DS_direction{t}{th})
                scatter(ROIx(select_id(c)), ROIy(select_id(c)), 'MarkerEdgeColor', colorRange(DS_direction{t}{th}(c), :), ...
                    'MarkerFaceColor', colorRange(DS_direction{t}{th}(c), :)); hold on
            end
            colormap hsv; colorbar('Ticks',0:0.125:1, 'TickLabels',{'0','45','90','135','180','225','270','315','360'}); 
            xlim([0 512]); ylim([0 512])
            title(['th', num2str(th), ' no fitting'])
            
            colorRange = hsv(361);
            subplot(1, 2, 2)
            select_id = id2{t}(fit_DS_id{t}{th});
            for c = 1 : length(fit_DS_direction{t}{th})
                scatter(ROIx(select_id(c)), ROIy(select_id(c)), 'MarkerEdgeColor', colorRange(1+uint16(fit_DS_direction{t}{th}(c)), :), ...
                    'MarkerFaceColor', colorRange(1+uint16(fit_DS_direction{t}{th}(c)), :)); hold on
            end
            title(['th', num2str(th), ' bimodal fitting'])
            colormap hsv; colorbar('Ticks',0:0.125:1, 'TickLabels',{'0','45','90','135','180','225','270','315','360'}); 
            xlim([0 512]); ylim([0 512])
            saveas(h, [save_folder, save_tag, 'spatial_p', num2str(t), '_th', num2str(th), '.png'])
            
            
            % plot averaged amplitude across all cells with DSI > threshold
            clear avg_amp
            if length(DS_direction{t}{th}) > 5
                h = figure; set(h, 'position', [100 100 1100 500])
                subplot(1, 3, 1)
                select_id = id1{t}(DS_id{t}{th});           
                for c = 1 : length(DS_direction{t}{th})
                    amp(c, :) = respMatrix{t}(select_id(c), :);
                    amp(c, :) = amp(c, :) / max(respMatrix{t}(select_id(c), :));
                    polarplot(((0 : 1/4 : 2) - 1/2) * pi, [amp(c, :), amp(c, 1)], 'linewidth', 1, 'color', [0.7 0.7 0.7]); hold on
                end
                avg_amp = mean(amp);
                polarplot(((0 : 1/4 : 2) - 1/2) * pi, [avg_amp, avg_amp(1)], 'linewidth', 4, 'color', 'k'); 
                title(['th', num2str(th), ' no fitting'])

                subplot(1, 3, 2)
                select_id = id2{t}(fit_DS_id{t}{th});           
                for c = 1 : length(fit_DS_direction{t}{th})
                    amp(c, :) = respMatrix{t}(select_id(c), :);
                    amp(c, :) = amp(c, :) / max(respMatrix{t}(select_id(c), :));
                    polarplot(((0 : 1/4 : 2) - 1/2) * pi, [amp(c, :), amp(c, 1)], 'linewidth', 1, 'color', [0.7 0.7 0.7]); hold on
                end
                avg_amp = mean(amp);
                polarplot(((0 : 1/4 : 2) - 1/2) * pi, [avg_amp, avg_amp(1)], 'linewidth', 4, 'color', 'k'); 
                title(['th', num2str(th), ' bimodal fitting'])

                subplot(1, 3, 3)
                select_id = id2{t};           
                for c = 1 : length(select_id)
                    amp(c, :) = respMatrix{t}(select_id(c), :);
                    amp(c, :) = amp(c, :) / max(respMatrix{t}(select_id(c), :));
                    polarplot(((0 : 1/4 : 2) - 1/2) * pi, [amp(c, :), amp(c, 1)], 'linewidth', 1, 'color', [0.7 0.7 0.7]); hold on
                end
                avg_amp = mean(amp);
                polarplot(((0 : 1/4 : 2) - 1/2) * pi, [avg_amp, avg_amp(1)], 'linewidth', 4, 'color', 'k'); 
                title('all responding cells')
                saveas(h, [save_folder, save_tag, 'avgAmp_p', num2str(t), '_th', num2str(th), '.png'])
            end
        end
        
      
        
        
        % plot the cells with 10% highest DSI        
        DSI_sort = sort(DSI(id1{t}, t));
        fit_DSI_sort = sort(DSI_fit(id2{t}, t));
        if length(DSI_sort) > 10
            sorted_id{t} = find(DSI(id1{t}, t) >= DSI_sort(uint16(length(DSI_sort)*0.9)));  
            fit_sorted_id{t} = find(DSI_fit(id2{t}, t) >= fit_DSI_sort(uint16(length(fit_DSI_sort)*0.9)));
            h = figure; set(h, 'position', [100 100 1100 500])
            subplot(1, 2, 1)
            DS_direction_sorted{t}{1} = directionId(id1{t}(sorted_id{t}), t);
            [counts_d, centers_d] = hist(DS_direction_sorted{t}{1}, 1:8);
            proportion_d = counts_d / length(DS_direction_sorted{t}{1});
            polarplot(((0 : 1/4 : 2) - 1/2) * pi, [proportion_d, proportion_d(1)], 'linewidth', 4)
            set(gca,'linew',2, 'FontSize', 15)
            saveas(h, [save_folder, save_tag, 'distribution_top10percent.png'])
            
            subplot(1, 2, 2)
            DS_direction_sorted{t}{2} = fit_results{t}(1, id2{t}(fit_sorted_id{t}));
            [counts_d, centers_d] = hist(DS_direction_sorted{t}{2}, 0:45:315);
            proportion_d = counts_d / length(DS_direction_sorted{t}{2});
            polarplot(((0 : 1/4 : 2) - 1/2) * pi, [proportion_d, proportion_d(1)], 'linewidth', 4)
            set(gca,'linew',2, 'FontSize', 15)
            saveas(h, [save_folder, save_tag, 'distribution_top10percent.png'])
        
        
            colorRange = hsv(8);
            h = figure; set(h, 'position', [100 100 1100 500])
            subplot(1, 2, 1)
            select_id = id1{t}(sorted_id{t});
            for c = 1 : length(DS_direction_sorted{t}{1})
                scatter(ROIx(select_id(c)), ROIy(select_id(c)), 'MarkerEdgeColor', colorRange(uint16(DS_direction_sorted{t}{1}(c)), :), ...
                    'MarkerFaceColor', colorRange(uint16(DS_direction_sorted{t}{1}(c)), :)); hold on
            end
            colormap jet; colorbar('Ticks',0:0.125:1, 'TickLabels',{'0','45','90','135','180','225','270','315','360'});
            xlim([0 512]); ylim([0 512])
            
            subplot(1, 2, 2)
            colorRange = hsv(361);
            select_id = id2{t}(fit_sorted_id{t});
            for c = 1 : length(DS_direction_sorted{t}{2})
                scatter(ROIx(select_id(c)), ROIy(select_id(c)), 'MarkerEdgeColor', colorRange(1+uint16(DS_direction_sorted{t}{2}(c)), :), ...
                    'MarkerFaceColor', colorRange(1+uint16(DS_direction_sorted{t}{2}(c)), :)); hold on
            end
            colormap jet; colorbar('Ticks',0:0.125:1, 'TickLabels',{'0','45','90','135','180','225','270','315','360'});
            xlim([0 512]); ylim([0 512])
            saveas(h, [save_folder, save_tag, 'spatial_top10percent.png'])
        end        
      
    end  
    
    
    
    % time to peak
    h = figure;
    subplot(2, 1, 1)
    hist(time2peak{t}, 10); hold on
    line([median(time2peak{t}), median(time2peak{t})], [0 5], 'color', 'r', 'linewidth', 2)
    title(['peak median=', num2str(median(time2peak{t}))])
    subplot(2, 1, 2)
    hist(time2slope{t}, 10); hold on
    line([median(time2slope{t}), median(time2slope{t})], [0 5], 'color', 'r', 'linewidth', 2)
    title(['slope median=', num2str(median(time2slope{t}))])    
    saveas(h, [save_folder, save_tag, 'respTime.png'])
    
    
    
    % tuning width
    h = figure;
    hist(tuning{t}, 10); hold on
    line([median(tuning{t}), median(tuning{t})], [0 5], 'color', 'r', 'linewidth', 2)
    title(['peak median=', num2str(median(tuning{t}))])
    saveas(h, [save_folder, save_tag, 'tuningWidth.png'])
end



h = figure; set(h, 'position', [0 0 1000 500])
subplot(1, 3, 1)
plot(1:3, resp_Percent, 'linewidth', 2); xlim([0 4]); ylim([0 1]); title('response percentage')

subplot(1, 3, 2)
plot(1:3, DS_percent(1, :), 'color', 'r', 'linewidth', 2); hold on
plot(1:3, DS_percent(2, :), 'color', 'k', 'linewidth', 2);
legend({'th1', 'th2'})
xlim([0 4]); ylim([0 1]); title('DS percentage')

subplot(1, 3, 3)
plot(OS_percent(1, :), 'color', 'r', 'linewidth', 2); hold on
plot(OS_percent(2, :), 'color', 'k', 'linewidth', 2);
legend({'th1', 'th2'})
xlim([0 4]); ylim([0 1]); title('OS percentage')
saveas(h, [save_folder, save_tag, 'response_percentage.png'])
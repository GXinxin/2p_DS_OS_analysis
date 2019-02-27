clear; clc;

tag = '3s_global';
% save_path = 'E:\Lab\Data\2p\summaryData\saline\';
% data_path = 'E:\Lab\Data\2p\saline';
save_path = 'E:\Lab\Data\2p\saline\kz_slc17a_ai162_P15_saline_2p_02132019_1\';
data_path = 'E:\Lab\Data\2p\saline\kz_slc17a_ai162_P15_saline_2p_02132019_1\';

cd(data_path)



% trial info, detection window
stimDuration = 10; % 10s grating
stimStatic = 5; % 5s static grating
afterStaticStim = 5; % 5s static grating after drifting grating

detectWd = 3; % use the first 3s to detect drift grating response
sampleFreq = 15;

timeStep1 = 16; % trace start 15 frames before onset of each trial
timeStep2 = timeStep1 + stimStatic * sampleFreq;
timeStep3 = timeStep2 + detectWd * sampleFreq;
timeStep4 = timeStep2 + (stimDuration - stimStatic) * sampleFreq;
timeStep5 = timeStep4 + afterStaticStim * sampleFreq;

  

% a_list = dir(fullfile('*2p*'));
a_list = dir(fullfile('*DS*'));

for a_id = 3% : length(a_list)

    clear isResponse isResponse2 meanResponse DSI maxResponse meanResponse_static
    cd(a_list(a_id).name)
    save_tag = ['individual_', num2str(detectWd), 's_'];
    save_folder = [save_path, a_list(a_id).name, '\'];
    mkdir(save_folder)


    % load responseTrace 
    flist = dir(fullfile('responseTrace*.mat'));
    if length(flist) == 1
        load(flist(1).name)
    else
        for f = 1 : length(flist)
            load(flist(f).name)
            if f == 1
                responseTrace_total = responseTrace;
            else
                for angle = 1 : length(responseTrace)
                    responseTrace_total{angle} = cat(3, responseTrace_total{angle}, responseTrace{angle});
                end
            end     
        end
        responseTrace = responseTrace_total;
    end


    % check cell firing, if skewness < 0.5, mark the cell as
    % inactive
    skew_mat = skewness(dF);
    inactive_id = find(skew_mat < 0.5);
    percent_inactive = length(inactive_id) / length(skew_mat)
    active_id = find(skew_mat >= 0.5);

    new_dF = dF(:, active_id);
    for d = 1 : 8
        respTrace{d} = responseTrace{d}(:, active_id, :);
    end



    % convert directions if right eye stimulated, direction 3 and 7 is
    % the same between left and right eyes
    if ~isempty(strfind(a_list(a_id).name, 'RE'))
        rightEye = 1;
        tmp = respTrace;
        respTrace{1} = tmp{5};
        respTrace{2} = tmp{4};
        respTrace{4} = tmp{2};
        respTrace{5} = tmp{1};
        respTrace{6} = tmp{8};
        respTrace{8} = tmp{6};
    else
        rightEye = 0;
    end



    for d = 1 : 8

        maxStatic{d} = squeeze(max(respTrace{d}(timeStep1 : timeStep2 - 1, :, :)));
        maxDrifting{d} = squeeze(max(respTrace{d}(timeStep2 : timeStep3, :, :)));

        meanStatic{d} = squeeze(mean(respTrace{d}(timeStep1 : timeStep2 - 1, :, :)));
        maxDrifting2{d} = squeeze(max(respTrace{d}(timeStep2 : timeStep4, :, :)));
        meanDrifting2{d} = squeeze(mean(respTrace{d}(timeStep2 : timeStep4, :, :)));

        isResponse(:, d) = ttest(maxDrifting{d}', maxStatic{d}', 'Tail', 'right'); % matched 2-sample
        isResponse2(:, d) = ttest2(maxDrifting{d}', maxStatic{d}', 'Tail', 'right'); % 2-sample   

        meanResponse(:, d) = mean(meanDrifting2{d}, 2); % - meanStatic{d}, 2);
        maxResponse(:, d) = mean(maxDrifting2{d}, 2); % averaged maximum firing rate for each direction
        meanResponse_static(:, d) = mean(maxStatic{d}, 2); % - meanStatic{d}, 2);
    end

%         respondId2 = find(sum(isResponse, 2) > 0);
    respondId2 = find(sum(isResponse2, 2) > 0);


    % percent of cells responding to any drift gratings
    percent(1) = length(respondId2) / size(isResponse, 1);
    percent


    % compute the DSI of responding cells
    responseMatrix = meanResponse(respondId2, :);
    responseMatrix_max = maxResponse(respondId2, :);
    maxDS = max(responseMatrix_max, [], 2);
    responseMatrix_static = meanResponse_static(respondId2, :);
    maxOS_static = max(responseMatrix_static, [], 2);
    
    DSI = [];
    DSI_global = [];
    OSI_global = [];
    directionId = [];
    for ii = 1 : length(maxDS)
        directionId(ii) = find(responseMatrix_max(ii, :) == maxDS(ii));
        oppositeId(ii) = mod(directionId(ii) + 4, 8);
        orthoId(ii, 1) = mod(directionId(ii) + 2, 8);
        orthoId(ii, 2) = mod(directionId(ii) - 2, 8);
        if oppositeId(ii) == 0
            oppositeId(ii) = 8;
        end
        orthoId(ii, orthoId(ii, :) == 0) = 8;


        DSI(ii) = (responseMatrix_max(ii, directionId(ii)) - responseMatrix_max(ii, oppositeId(ii)))...
            / (responseMatrix_max(ii, directionId(ii)) + responseMatrix_max(ii, oppositeId(ii)));

        R_ortho = mean(responseMatrix(ii, orthoId(ii, :)));
        OSI(ii) = (responseMatrix_max(ii, directionId(ii)) - R_ortho)...
            / (responseMatrix_max(ii, directionId(ii)) + R_ortho);
        
        
        % use static grating to compute OS
        orientationId(ii) = find(responseMatrix_static(ii, :) == maxOS_static(ii));
        orthoId_static(ii, 1) = mod(orientationId(ii) + 2, 8);
        orthoId_static(ii, 2) = mod(orientationId(ii) - 2, 8);
        orthoId_static(ii, orthoId_static(ii, :) == 0) = 8;
        R_ortho_static = mean(responseMatrix_static(ii, orthoId_static(ii, :)));
        OSI_static(ii) = (responseMatrix_static(ii, orientationId(ii)) - R_ortho_static)...
            / (responseMatrix_static(ii, orientationId(ii)) + R_ortho_static);

        
        % global DSI
%         tmp = responseMatrix(ii, :) /  max(responseMatrix(ii, :));
        tmp = responseMatrix_max(ii, :) /  max(responseMatrix_max(ii, :));
        DSI_global(ii) = abs(sum(tmp .* exp(1i * pi/4 * (1:8))) / (sum(tmp)));
        OSI_global(ii) = abs(sum(tmp .* exp(2 * 1i * pi/4 * (1:8))) / (sum(tmp)));
        
        tmp_static = responseMatrix_static(ii, :) /  max(responseMatrix_static(ii, :));
        OSI_global_static(ii) = abs(sum(tmp_static .* exp(2 * 1i * pi/4 * (1:8))) / (sum(tmp_static)));
    end


    % percentage of cells having DSI larger than .3
    DS_id{1} = find(DSI > 0.3);
    percent(2) = length(DS_id{1}) / length(DSI);
    DS_direction{1} = directionId(DS_id{1});

    DS_id{2} = find(DSI > 0.5);
%     percent(2) = length(DS_id{2}) / length(DSI);
    DS_direction{2} = directionId(DS_id{2});

    median_DSI = median(DSI);

    save([save_folder, 'DSI_summary_window', num2str(detectWd), '.mat'], 'directionId', 'orientationId', 'percent', 'DS_direction', ...
        'median_DSI', 'DSI', 'OSI', 'respondId2', 'DS_id', 'DSI_global', 'OSI_global', 'percent_inactive', 'OSI_static', 'OSI_global_static')




    % plot DS distribution and polar plots for each animal
    colorVector = copper(4);

    h = figure;
    [counts, centers] = hist(DSI, 0.08:0.08:0.8);
    proportion = counts / length(DSI);
    bar(centers, proportion, 'barWidth', 0.9)
    line(median_DSI * [1 1], [0 max(proportion)], 'color', 'r', 'lineWidth', 2)
    xlim([0 0.9])
%         ylim([0 0.3])   
    ylabel('proportion', 'FontSize', 12, 'FontWeight', 'Bold')
    xlabel('DSI', 'FontSize', 12, 'FontWeight', 'Bold')
    box off
    saveas(h, [save_folder, save_tag, 'proportion1.png'])



    h = figure;
    [counts_d, centers_d] = hist(DS_direction{1}, 1:8);
    proportion_d = counts_d / length(DS_direction{1});
    polarplot(((0 : 1/4 : 2) - 1/2) * pi, [proportion_d, proportion_d(1)], 'linewidth', 4)
    rlim([0, 0.5])
    set(gca,'linew',2, 'FontSize', 22)
    saveas(h, [save_folder, save_tag, 'distribution.png'])


    cd ..
end







% added bimodal distribution fitting for drift grating
% changed criterion for cell responsiveness



clear; clc;
tag = 'all';

fileInfo = readtext('E:\Lab\Data\2p\analysisInputs\DG_summaryEachFOVInfo_8s.txt', ' ');
for f = 1 : size(fileInfo, 1)
    compute_individual_DS(fileInfo{f, 1}, fileInfo{f, 2}, fileInfo{f, 3}, fileInfo{f, 4});
end


% compute DSI, OSI, tuning width, etc...
function compute_individual_DS(data_path, data_group, data_tag, rightEye)

disp(data_path(22:end))
switch data_group
    case 7
        save_path = ['E:\Lab\Data\2p\tra2b\v1\summary', data_path(26:end)]; % h2b inj in V1
    case 14
        save_path = ['E:\Lab\Data\2p\tra2b\summary\gabazine', data_path(30:end)]; % tra2b p15 gabazine
    case 15
        save_path = ['E:\Lab\Data\2p\tra2b\summary\saline', data_path(30:end)]; % tra2b p15 gabazine
    case 16
        save_path = ['E:\Lab\Data\2p\tra2b\summary\p21\gabazine', data_path(30:end)]; % tra2b p15 gabazine
    case 17
        save_path = ['E:\Lab\Data\2p\tra2b\summary\p21\saline', data_path(30:end)]; % tra2b p15 gabazine
end


% if there are multiple recording session saved under one single folder,
% then separate these folders
if strcmp(data_tag, 'acq*')
    slashId = strfind(data_path, '\');
    save_path = [save_path, '\', data_path(slashId(end)+1 : end)];
end


cd(data_path)



% trial info, detection window
stimDuration = 16; % 10s grating
stimStatic = 8; % 5s static grating
afterStaticStim = 8; % 5s static grating after drifting grating

detectWd = 10; % use the first 3s to detect drift grating response
sampleFreq = 15;

timeStep(1) = 16; % trace start 15 frames before onset of each trial
timeStep(2) = timeStep(1) + stimStatic * sampleFreq; % offset of static gratings
timeStep(3) = timeStep(2) + detectWd * sampleFreq; % offset of drifting gratings
% try 2 + timeStep(3) in 8s condition due to the slow dynamics of H2b
timeStep(4) = timeStep(2) + stimDuration * sampleFreq; % offset of drift gratings
% timeStep(5) = timeStep(4) + (afterStaticStim-2) * sampleFreq; % end of the trial


% load data from each subfolder
a_list = dir(fullfile(data_tag));

if length(a_list) >= 1
    for a_id = 1 : length(a_list)
        a_id
        [data_path, a_list(a_id).name]
        clear isResponse isResponse2 meanResponse DSI maxResponse meanResponse_static meanAmp preActive activated excludeActivation...
            is_respond0 excludeActiveCell is_respond respMatrix maxDS maxOS excludeActivation maxFrame is_selective anova_p...
            fit_rsquare fit_error iscell active_id
        cd([data_path, a_list(a_id).name])
        save_tag = ['tra2b_', num2str(detectWd), 's_zscore_'];
        save_folder = [save_path, a_list(a_id).name, '\'];
        mkdir(save_folder)
        
        
        
        % load responseTrace, concatenate if had multiple sessions
        flist = dir(fullfile('acq*'));
        
        
        for ff = 1 : length(flist)
            cd(flist(ff).name)
            load('responseTrace_noZscore.mat')
            if ff == 1
                responseTrace_total = responseTrace;
                iscell_total = add_iscell;
            else
                for angle = 1 : length(responseTrace)
                    responseTrace_total{angle} = cat(3, responseTrace_total{angle}, responseTrace{angle});
                    iscell_total = iscell_total .* add_iscell;
                end
            end
            cd ..
        end
        
        iscell_final = iscell_total(is_cell > 0);
        responseTrace = responseTrace_total;
        for angle = 1 : length(responseTrace)
            responseTrace{angle} = responseTrace{angle}(:, iscell_final>0, :);
        end
        dF = dF(:, iscell_final>0);
        
        
        cd ..
        
        
        
        % check cell firing, if skewness < 0.5, mark the cell as
        % inactive
        skew_mat = skewness(dF);
        inactive_id = find(skew_mat < 0.5);
        percent_inactive = length(inactive_id) / length(skew_mat);
        active_id = find(skew_mat >= 0.5);
        
        new_dF = dF(:, active_id);
        for d = 1 : 8
            respTrace{d} = responseTrace{d}(:, active_id, :);
        end
        no_cell = size(new_dF, 2);
        
        
        
        % compute mean and std for each cell
        for c = 1 : no_cell
            std_cell = std(new_dF);
            mean_cell = mean(new_dF);
        end
        
        
        
        % convert directions if right eye stimulated, direction 3 and 7 is
        % the same between left and right eyes
        if rightEye
            tmp = respTrace;
            respTrace{1} = tmp{5};
            respTrace{2} = tmp{4};
            respTrace{4} = tmp{2};
            respTrace{5} = tmp{1};
            respTrace{6} = tmp{8};
            respTrace{8} = tmp{6};
        end
        
        
        
        % compute averaged response across all directions for each cell, and
        % find the frame where the peak response occurred
        for d = 1 : 8
            for t = 1 : 3
                meanTraces{t}{d} = squeeze(mean(respTrace{d}(timeStep(t) : timeStep(t+1) - 1, :, :), 3));
                meanResponse(:, d, t) = max(meanTraces{t}{d}); % find the peak of the averaged response during each stimuli period
                
                for c = 1 : size(meanResponse, 1) % find the frame index where the peak of averaged response occurred
                    maxFrame(c, d, t) = find(meanTraces{t}{d}(:, c) == meanResponse(c, d, t), 1);
                end
                preActive(:, d, t) = maxFrame(:, d, t) < 10; % label cells with activation occurred prior each stimuli period: whether peak is earlier than 10 frames
                meanAmp{t}{d} = squeeze(mean(respTrace{d}(timeStep(t) : timeStep(t+1) - 1, :, :)));
                %                 meanAmp{t}(:, :, d) = squeeze(mean(respTrace{d}(timeStep(t) : timeStep(t+1) - 1, :, :)));
            end
        end
        
        
        
        % check whether cells respond to each stimuli period: static, drifting,
        % static after drifting: if peak value of any direction > 1 std + mean, then responsive
        for t = 1 : 3
            activated(:, :, t) = meanResponse(:, :, t) > repmat(std_cell' + mean_cell', 1, 8);
            is_respond0(:, t) = sum(activated(:, :, t), 2) > 0;
            
            excludeActivation(:, :, t) = preActive(:, :, t) .* activated(:, :, t); % exclude activations that occur before each stimuli period
            excludeActiveCell(:, t) = (sum(activated(:, :, t), 2) - sum(excludeActivation(:, :, t), 2)) == 0;
            is_respond(:, t) = is_respond0(:, t) .* (1 - excludeActiveCell(:, t));
            
            resp_Percent(t) = sum(is_respond(:, t)) / length(is_respond(:, t)); % percentage of cells that respond to each stimuli
        end
        
        
        
        % compute the DSI of responding cells for each stimuli period
        % separately
        DSI = [];
        DSI_global = [];
        OSI_global = [];
        directionId = [];
        
        for t = 1 : 3
            for d = 1 : 8
                respMatrix{t}(:, d) = double(squeeze(mean(meanAmp{t}{d}, 2)));
            end
            respMatrix_axis{t} = (respMatrix{t}(:, 1:4) + respMatrix{t}(:, 5:8)) / 2;
            maxDS(:, t) = max(respMatrix{t}, [], 2);
            maxOS(:, t) = max(respMatrix_axis{t}, [], 2);
            
            
            for ii = 1 : no_cell
                
                % test whether anisotropy by one-way ANOVA
                tmp = [];
                group = [];
                for d = 1 : 8
                    tmp = [tmp, meanAmp{t}{d}(ii, :)];
                    group = [group, repmat({num2str(d)}, 1, length(meanAmp{t}{d}(ii, :)))];
                end
                anova_p(ii, t) = anova1(tmp, group, 'off');
                
                directionId(ii, t) = find(respMatrix{t}(ii, :) == maxDS(ii, t)); % note the DSI during static grating is meaningless, just include for coding convenience
                oppositeId(ii, t) = mod(directionId(ii, t) + 4, 8);
                orientationId(ii, t) = find(respMatrix_axis{t}(ii, :) == maxOS(ii, t));
                orthoId(ii, t) = mod(orientationId(ii, t) + 2, 4);
                if oppositeId(ii, t) == 0
                    oppositeId(ii, t) = 8;
                end
                if orthoId(ii, t) == 0
                    orthoId(ii, t) = 4;
                end
                
                
                % compute DSI and OSI
                DSI(ii, t) = (respMatrix{t}(ii, directionId(ii, t)) - respMatrix{t}(ii, oppositeId(ii, t)))...
                    / (respMatrix{t}(ii, directionId(ii, t)) + respMatrix{t}(ii, oppositeId(ii, t)));
                
                OSI(ii, t) = (respMatrix_axis{t}(ii, orientationId(ii, t)) - respMatrix_axis{t}(ii, orthoId(ii, t)))...
                    / (respMatrix_axis{t}(ii, orientationId(ii, t)) + respMatrix_axis{t}(ii, orthoId(ii, t)));
                
                
                % global DSI and OSI
                tmp = respMatrix{t}(ii, :) /  max(respMatrix{t}(ii, :));
                DSI_global(ii, t) = abs(sum(tmp .* exp(1i * pi/4 * (1:8))) / (sum(tmp)));
                OSI_global(ii, t) = abs(sum(tmp .* exp(2 * 1i * pi/4 * (1:8))) / (sum(tmp)));
                
                
                % fit bimodal gaussian for computing orientation selectivity
                [fitCoeff{ii, t}, gof{ii, t}] = fitBimodalGaussian(respMatrix{t}(ii, :), 0);
                fit_results{t}(1, ii) = fitCoeff{ii, t}.pref;
                fit_results{t}(2, ii) = fitCoeff{ii, t}.rp;
                fit_results{t}(3, ii) = fitCoeff{ii, t}.ro;
                fit_results{t}(4, ii) = fitCoeff{ii, t}.rc;
                fit_results{t}(5, ii) = 2 * sqrt(2*log(2)) * fitCoeff{ii, t}.sigma; % half the width of the gaussian at half its height
                
                
                % apply the fitted gaussian to compute DSI and OSI
                rc = fit_results{t}(4, ii);
                rp = fit_results{t}(2, ii);
                ro = fit_results{t}(3, ii);
                pref = fit_results{t}(1, ii);
                sigma = fitCoeff{ii, t}.sigma;
                gaussEqn = @(x) rc + rp * exp(-(min([abs(x-pref), abs(360-x+pref), abs(x-pref+360)], [], 2)/(sqrt(2)*sigma))^2) + ro * exp(-(min([abs(x-pref+180), abs(180-x+pref), abs(x-pref+540)], [], 2)/(sqrt(2)*sigma))^2)';
                DSI_fit(ii, t) = (gaussEqn(pref) - gaussEqn(pref+180)) / (gaussEqn(pref) + gaussEqn(pref+180)); % DSI computed by fitted curves
                ori_pref = (gaussEqn(pref) + gaussEqn(pref+180));
                ori_ortho = (gaussEqn(pref+90) + gaussEqn(pref-90));
                OSI_fit(ii, t) = (ori_pref - ori_ortho) / (ori_pref + ori_ortho);
                
                % error and r-square of fit
                fit_error(ii, t) = gof{ii, t}.sse;
                fit_rsquare(ii, t) = gof{ii, t}.rsquare;
            end
        end
        % good fits
        is_fit = (fit_error < 0.4) .* (fit_rsquare > 0.6);
        % is anisotropy by one-way ANOVA
        is_selective = anova_p < 0.05;
        
        
        
        % summarize DSI and OSI for each stimuli period
        for t = 1 : 3
            % percentage of cells having DSI larger than .3 or .5
            thresh = [0.3 0.5];
            for th = 1 : 2
                % measured results, include cells with peak of averaged trace
                % larger than 1 std + mean, and has significant p-value from
                % one-way anova
                id1{t} = find(is_respond(:, t) .* is_selective(:, t) > 0);
                DS_id{t}{th} = find(DSI(id1{t}, t) > thresh(th));
                DS_percent(th, t) = length(DS_id{t}{th}) / length(id1{t});
                DS_direction{t}{th} = directionId(id1{t}(DS_id{t}{th}), t); % was directionId(DS_id{t}{th}, t)
                
                OS_id{t}{th} = find(OSI(id1{t}, t) > thresh(th));
                OS_percent(th, t) = length(OS_id{t}{th}) / length(id1{t});
                OS_direction{t}{th} = orientationId(id1{t}(OS_id{t}{th}), t); % was directionId(OS_id{t}{th}, t)
                
                
                % fitted results, include cells with peak of averaged trace
                % larger than 1 std + mean, has good bimodal gaussian fit, and
                % has p<0.05 with one-way ANOVA
                id2{t} = find(is_respond(:, t) .* is_fit(:, t) .* is_selective(:, t) > 0);
                fit_DS_id{t}{th} = find(DSI_fit(id2{t}, t) > thresh(th));
                fit_DS_percent(th, t) = length(fit_DS_id{t}{th}) / length(id2{t});
                fit_DS_direction{t}{th} = fit_results{t}(1, id2{t}(fit_DS_id{t}{th})); % was fit_results{t}(1, fit_DS_id{t}{th})
                tmp = find(fit_DS_direction{t}{th} <= 0);
                fit_DS_direction{t}{th}(tmp) = fit_DS_direction{t}{th}(tmp) + 360;
                
                fit_OS_id{t}{th} = find(OSI_fit(id2{t}, t) > thresh(th));
                fit_OS_percent(th, t) = length(fit_OS_id{t}{th}) / length(id2{t});
                fit_OS_direction{t}{th} = fit_results{t}(1, id2{t}(fit_OS_id{t}{th})); % was fit_results{t}(1, fit_OS_id{t}{th})
                tmp = find(fit_OS_direction{t}{th} > 180);
                fit_OS_direction{t}{th}(tmp) = fit_OS_direction{t}{th}(tmp) - 180;
                tmp = find(fit_OS_direction{t}{th} <= 0);
                fit_OS_direction{t}{th}(tmp) = fit_OS_direction{t}{th}(tmp) + 180;
            end
            
            median_values(t, 1) = median(DSI(id1{t}, t));
            median_values(t, 2) = median(DSI_fit(id2{t}, t));
            median_values(t, 3) = median(DSI_global(id1{t}, t));
            median_values(t, 4) = median(OSI(id1{t}, t));
            median_values(t, 5) = median(OSI_fit(id2{t}, t));
            median_values(t, 6) = median(OSI_global(id1{t}, t));
            tuning{t} = fit_results{t}(5, id2{t});
            
            
            
            % compute the resp onset time during drift grating: time to slope and
            % time to peak
            count = 1;
            time2peak{t} = [];
            time2slope{t} = [];
            for c = 1 : length(DS_id{t}{1})
                trace = meanTraces{t}{directionId(DS_id{t}{1}(c), t)}(:, DS_id{t}{1}(c));
                diffTrace = trace(3:end) - trace(1:end-2);
                time2peak{t}(count) = find(trace == max(trace));
                time2slope{t}(count) = find(diffTrace == max(diffTrace)) + 1;
                count = count + 1;
            end
            
        end
        
        
        % save data for each FOV
        save([save_folder, 'DSI_summary_window', num2str(detectWd), '.mat'], 'directionId', 'orientationId', 'DS_direction', 'anova_p', ...
            'median_values', 'DSI', 'OSI', 'DSI_fit', 'OSI_fit', 'DS_id', 'DS_percent', 'DS_direction', 'OS_id', 'OS_percent', 'OS_direction', 'fit_DS_id', 'fit_DS_percent', 'fit_DS_direction', ...
            'fit_OS_id', 'fit_OS_percent', 'fit_OS_direction', 'DSI_global', 'OSI_global', 'percent_inactive', 'DSI_fit', 'OSI_fit', 'is_respond', 'id1', 'id2',...
            'is_fit', 'is_selective', 'time2peak', 'time2slope', 'respMatrix', 'respMatrix_axis', 'tuning', 'fitCoeff', 'gof', 'fit_results')
        
        
        
        
        % Plots
        cd([data_path, a_list(a_id).name])
        plot_individualResults
        cd ..
        close all
    end
end
end








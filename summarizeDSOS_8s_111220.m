

% path = 'E:\Lab\Data\2p\summaryData\newDataSet_042219\deep\';
% path = 'E:\Lab\Data\2p\summaryData\newDataSet_042219\deep\sameDeeper\';
% path = 'E:\Lab\Data\2p\tra2b\summary\';
% path = 'E:\Lab\Data\2p\tra2b\summary\onlyCenterData\';
% path = 'E:\Lab\Data\2p\tra2b\summary\onlyCenterData\addCorrData\';
path = 'E:\Lab\Data\wavePaperRevision\acuteGabazine\tra2b_p15_acute\';
% path = 'E:\Lab\Data\wavePaperRevision\acuteGabazine\tra2b_frmd7_p15\';
% path = 'E:\Lab\Data\wavePaperRevision\acuteGabazine\tra2b_p15_acute_noninj\';

cd(path)

% tag = {'saline', 'gabazine_20x', 'gabazine_100x'};
% tag = {'saline', 'gabazine', 'acute', 'exp'};
% tag = {'ctl', 'exp'};
% tag = {'saline', 'gabazine'};
% tag = {'noninj', 'acute'}
tag = {'saline', 'acute'};
% save_tag = 'p15_tra2b_';


thresh = [0.3 0.5]; % DSI th = 0.3 or th = 0.5
include_th = [0.2 0.3];
is_smooth = 0;
% colorVector = copper(4);
% colorVector = [180,160,0; 105,105,105] / 255; % saline vs gabazine
% colorVector = [180 200 50; 171 149 132] / 255; % control vs frmd7
% colorVector = [180 200 50; 30 10 130] / 255; % non-inj vs acute gabazine
% colorVector = [180,160,0; 105,105,105; 30 10 130; 171 149 132] / 255; % saline, chronic gabazine, acute gabazine, frmd7
colorVector = [180,160,0; 30 10 130] / 255; % saline, acute gabazine
colorVector_line = [255, 255, 255; 105,105,105] / 255; 
% colorVector = [204,205,42; 105,105,105]

% pixel_res = 410 / 512; % 410um for 512 pixels
pixel_res = 1; % new version corrected pixel to distance in previous steps
bin_um = 50; % bin size in um
bin_pixel = bin_um / pixel_res;


for ii = 1 %: length(include_th)
    if is_smooth
        fd_name = [num2str(ii), '_smooth'];
        mkdir(fd_name)
    else
        fd_name = [num2str(ii), '_noSmooth'];
        mkdir(fd_name)
    end
    save_tag = [fd_name, '\p15_tra2b_'];
    
    % initialization
%     tt = 3;
%     DSI_total = cell(tt, length(tag));
%     OSI_total = cell(tt, length(tag));
%     DSI_fit_total = cell(tt, length(tag));
%     OSI_fit_total = cell(tt, length(tag));
%     DSI_global_total = cell(tt, length(tag));
%     OSI_global_total = cell(tt, length(tag));
%     time2peak_total = cell(tt, length(tag));
%     time2slope_total = cell(tt, length(tag));
%     cell_no = cell(tt, length(tag));
%     resp_Percent_all = cell(tt, length(tag));
%     select_Percent_all = cell(tt, length(tag));
%     percent_DS = cell(tt, length(tag));
%     percent_OS = cell(tt, length(tag));
%     percent_DS_fit = cell(tt, length(tag));
%     percent_OS_fit = cell(tt, length(tag));
%     DS_direction_total = cell(tt, length(thresh), length(tag));
%     DS_direction_fit_total = cell(tt, length(thresh), length(tag));
%     OS_direction_total = cell(tt, length(thresh), length(tag));
%     OS_direction_fit_total = cell(tt, length(thresh), length(tag));
%     OS_DS_direction_total = cell(tt, length(thresh), length(tag));
%     OS_DS_direction_direction_total = cell(tt, length(thresh), length(tag));
%     OS_only_direction_total = cell(tt, length(thresh), length(tag));
%     fit_OS_DS_direction_total = cell(tt, length(thresh), length(tag));
%     fit_OS_only_direction_total = cell(tt, length(thresh), length(tag));
%     fit_OS_DS_direction_direction_total = cell(tt, length(thresh), length(tag));
%     fov_circ_var = cell(tt, length(thresh), length(tag));
%     total_pdist = cell(tt, length(thresh), length(tag));
%     individual_pdist = cell(tt, length(thresh), length(tag));
    %     aDSI_total = cell(3, 2);
    %     OSI_total = cell(3, 2);
    %     DSI_fit_total = cell(3, 2);
    %     OSI_fit_total = cell(3, 2);
    %     DSI_global_total = cell(3, 2);
    %     OSI_global_total = cell(3, 2);
    %     time2peak_total = cell(3, 2);
    %     time2slope_total = cell(3, 2);
    %     cell_no = cell(3, 2);
    %     resp_Percent_all = cell(3, 2);
    %     select_Percent_all = cell(3, 2);
    %     percent_DS = cell(3, 2);
    %     percent_OS = cell(3, 2);
    %     percent_DS_fit = cell(3, 2);
    %     percent_OS_fit = cell(3, 2);
    %     DS_direction_total = cell(3, 2, 2);
    %     DS_direction_fit_total = cell(3, 2, 2);
    %     OS_direction_total = cell(3, 2, 2);
    %     OS_direction_fit_total = cell(3, 2, 2);
    %     OS_DS_direction_total = cell(3, 2, 2);
    %     OS_DS_direction_direction_total = cell(3, 2, 2);
    %     OS_only_direction_total = cell(3, 2, 2);
    %     fit_OS_DS_direction_total = cell(3, 2, 2);
    %     fit_OS_only_direction_total = cell(3, 2, 2);
    %     fit_OS_DS_direction_direction_total = cell(3, 2, 2);
    %     fov_circ_var = cell(3, 2, 2);
    %     total_pdist = cell(3, 2, 2);
    %     individual_pdist = cell(3, 2, 2);
    
    
    for g = 2 : length(tag)
        
        cd([path, tag{g}])
        fdlist = dir(fullfile('*tra2b*'));
        
        count = 1
        
        is_include{g} = zeros(length(fdlist), 2);
        
        % for each animal
        for ff = 1 : length(fdlist)
            cd(fdlist(ff).name)
            sub_fdlist = dir(fullfile('*tra2b*'));
            
            % for each FOV per animal
            for sub_f = 1 : length(sub_fdlist)
                indexingM(:, count) = [ff; sub_f];
                cd(sub_fdlist(sub_f).name)
                if is_smooth
                    load('filtered_DSI_summary_window10.mat')
                else
                    load('DSI_summary_window10addCorr.mat')
                    %                     load('DSI_summary_window10useWeighted.mat')
                end
                
                % exclude FOVs with less than 20% responsive cells during
                % the drift grating period
                is_include{g}(ff, sub_f) = (sum(is_respond(:, 2)) / length(is_respond(:, 2))) > include_th(ii);
                if is_include{g}(ff, sub_f)
                    for t = 1 : 3
                        cell_no{t, g}(:, count) = [sum(is_respond(:, t)); length(is_respond(:, t)); length(id1{t})];
                        resp_Percent_all{t, g}(count) = sum(is_respond(:, t)) / length(is_respond(:, t)); % percentage of cells that respond to each stimuli
                        select_Percent_all{t, g}(count) = sum(is_respond(:, t) .* is_selective(:, t)) / length(is_respond(:, t));
                        percent_DS{t, g}(:, count) = DS_percent(:, t);
                        percent_OS{t, g}(:, count) = OS_percent(:, t);
                        percent_DS_fit{t, g}(:, count) = fit_DS_percent(:, t);
                        percent_OS_fit{t, g}(:, count) = fit_OS_percent(:, t);
                        percent_OS_DS{t, g}(:, count) = OS_DS_percent(:, t);
                        percent_OS_only{t, g}(:, count) = OS_only_percent(:, t);
                        percent_OS_DS_fit{t, g}(:, count) = fit_OS_DS_percent(:, t);
                        percent_OS_only_fit{t, g}(:, count) = fit_OS_only_percent(:, t);
                        
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
                            OS_DS_direction_total{t, th, g} = [OS_DS_direction_total{t, th, g}; OS_DS_direction{t}{th}];
                            OS_DS_direction_direction_total{t, th, g} = [OS_DS_direction_direction_total{t, th, g}; OS_DS_direction_direction{t}{th}];
                            OS_only_direction_total{t, th, g} = [OS_only_direction_total{t, th, g}; OS_only_direction{t}{th}];
                            fit_OS_DS_direction_total{t, th, g} = [fit_OS_DS_direction_total{t, th, g}, fit_OS_DS_direction{t}{th}];
                            fit_OS_only_direction_total{t, th, g} = [fit_OS_only_direction_total{t, th, g}, fit_OS_only_direction{t}{th}];
                            fit_OS_DS_direction_direction_total{t, th, g} = [fit_OS_DS_direction_direction_total{t, th, g}, fit_OS_DS_direction_direction{t}{th}];
                            
                            
                            
                            % also saves DS and OS properties for individual FOVs
                            DS_direction_individual{t, th, g}{count} = DS_direction{t}{th};
                            DS_direction_individual_fit{t, th, g}{count} = fit_DS_direction{t}{th};
                            OS_direction_individual{t, th, g}{count} = OS_direction{t}{th};
                            OS_direction_individual_fit{t, th, g}{count} = fit_OS_direction{t}{th};
                            OS_DS_direction_individual{t, th, g}{count} = OS_DS_direction{t}{th};
                            OS_DS_direction_direction_total_individual{t, th, g}{count} = OS_DS_direction_direction{t}{th};
                            OS_only_direction_total_individual{t, th, g}{count} = OS_only_direction{t}{th};
                            OS_DS_direction_total_individual_fit{t, th, g}{count} = fit_OS_DS_direction{t}{th};
                            OS_only_direction_total_individual_fit{t, th, g}{count} = fit_OS_only_direction{t}{th};
                            OS_DS_direction_direction_total_individual_fit{t, th, g}{count} = fit_OS_DS_direction_direction{t}{th};
                            
                            
                            % compute circular variance and percentage of
                            % directions/orientations along the opposite
                            % direction
                            for case_i = 1 : 11
                                switch case_i
                                    case 1
                                        tmp = OS_direction{t}{th} * pi/4;
                                        tmp2 = OS_direction{t}{th};
                                    case 2
                                        tmp = (fit_OS_direction{t}{th} * pi/180)';
                                        tmp2 = fit_OS_direction{t}{th};
                                    case 3
                                        tmp = DS_direction{t}{th} * pi/4;
                                        tmp2 = DS_direction{t}{th};
                                    case 4
                                        tmp = (fit_DS_direction{t}{th} * pi/180)';
                                        tmp2 = fit_DS_direction{t}{th};
                                    case 5
                                        tmp = OS_DS_direction{t}{th} * pi/4;
                                        tmp2 = OS_DS_direction{t}{th};
                                    case 6
                                        tmp = OS_only_direction{t}{th} * pi/4;
                                        tmp2 = OS_only_direction{t}{th};
                                    case 7
                                        tmp = OS_DS_direction_direction{t}{th} * pi/4;
                                        tmp2 = OS_DS_direction_direction{t}{th};
                                    case 8
                                        tmp = (fit_OS_DS_direction{t}{th} * pi/180)';
                                        tmp2 = fit_OS_DS_direction{t}{th};
                                    case 9
                                        tmp = (fit_OS_only_direction{t}{th} * pi/180)';
                                        tmp2 = fit_OS_only_direction{t}{th};
                                    case 10
                                        tmp = (fit_OS_DS_direction_direction{t}{th} * pi/180)';
                                        tmp2 = fit_OS_DS_direction_direction{t}{th};
                                    case 11
                                        tmp = OS_DS_direction_direction{t}{th} * pi/4;
                                        tmp2 = OS_DS_direction_direction{t}{th};
                                end
                                
                                if ~isempty(tmp)
                                    fov_circ_var{t, th, g}(count, case_i) = circ_var(tmp);
                                end
                                
                                
                                
                                if ~isempty(tmp)
                                    if ismember(case_i, [1 5 6 11])
                                        within_n = sum(tmp2 == 1);
                                    elseif ismember(case_i, [4 10])
                                        within_n = sum((tmp2 < 202.5) & (tmp2 > 157.5));
                                    elseif ismember(case_i, [2 8 9])
                                        within_n = sum((tmp2 < 22.5) + (tmp2 > 157.5));
                                    else
                                        within_n = sum(tmp2 == 5);
                                    end
                                    angle_proportion{t, th, g}(count, case_i) = within_n / length(tmp);
                                end
                                
                                
                                % store data for animal-based analysis
                                if t == 2 && th == 1
                                    switch case_i
                                        case 1
                                            tmp_animal{g}{case_i}{ff, sub_f} = OS_direction{t}{th} * pi/4;
                                            tmp_animal2{g}{case_i}{ff, sub_f} = OS_direction{t}{th};
                                        case 2
                                            tmp_animal{g}{case_i}{ff, sub_f} = (fit_OS_direction{t}{th} * pi/180)';
                                            tmp_animal2{g}{case_i}{ff, sub_f} = fit_OS_direction{t}{th};
                                        case 3
                                            tmp_animal{g}{case_i}{ff, sub_f} = DS_direction{t}{th} * pi/4;
                                            tmp_animal2{g}{case_i}{ff, sub_f} = DS_direction{t}{th};
                                        case 4
                                            tmp_animal{g}{case_i}{ff, sub_f} = (fit_DS_direction{t}{th} * pi/180)';
                                            tmp_animal2{g}{case_i}{ff, sub_f} = fit_DS_direction{t}{th};
                                        case 5
                                            tmp_animal{g}{case_i}{ff, sub_f} = OS_DS_direction{t}{th} * pi/4;
                                            tmp_animal2{g}{case_i}{ff, sub_f} = OS_DS_direction{t}{th};
                                        case 6
                                            tmp_animal{g}{case_i}{ff, sub_f} = OS_only_direction{t}{th} * pi/4;
                                            tmp_animal2{g}{case_i}{ff, sub_f} = OS_only_direction{t}{th};
                                        case 7
                                            tmp_animal{g}{case_i}{ff, sub_f} = OS_DS_direction_direction{t}{th} * pi/4;
                                            tmp_animal2{g}{case_i}{ff, sub_f} = OS_DS_direction_direction{t}{th};
                                        case 8
                                            tmp_animal{g}{case_i}{ff, sub_f} = (fit_OS_DS_direction{t}{th} * pi/180)';
                                            tmp_animal2{g}{case_i}{ff, sub_f} = fit_OS_DS_direction{t}{th};
                                        case 9
                                            tmp_animal{g}{case_i}{ff, sub_f} = (fit_OS_only_direction{t}{th} * pi/180)';
                                            tmp_animal2{g}{case_i}{ff, sub_f} = fit_OS_only_direction{t}{th};
                                        case 10
                                            tmp_animal{g}{case_i}{ff, sub_f} = (fit_OS_DS_direction_direction{t}{th} * pi/180)';
                                            tmp_animal2{g}{case_i}{ff, sub_f} = fit_OS_DS_direction_direction{t}{th};
                                        case 11
                                            tmp_animal{g}{case_i}{ff, sub_f} = (OS_DS_direction_direction{t}{th} * pi/180)';
                                            tmp_animal2{g}{case_i}{ff, sub_f} = OS_DS_direction_direction{t}{th};
                                    end
                                    
                                    selected_directionId = directionId(id1{2}, 2);
                                    for d = 1 : 8
                                        cell_id = find(selected_directionId == d);
                                        tmp = DSI(id1{2}, 2);
                                        DSI_perDir{g}{d}{count} = tmp(cell_id);
                                        tmp = DSI_fit(id1{2}, 2);
                                        DSI_fit_perDir{g}{d}{count} = tmp(cell_id);
                                        tmp = DSI_global(id1{2}, 2);
                                        DSI_global_perDir{g}{d}{count} = tmp(cell_id);
                                    end
                                    
                                    selected_orientationId = orientationId(id1{2}, 2);
                                    for o = 1 : 4
                                        cell_id = find(selected_orientationId == o);
                                        tmp = OSI(id1{2}, 2);
                                        OSI_perOri{g}{o}{count} = tmp(cell_id);
                                        tmp = OSI_fit(id1{2}, 2);
                                        OSI_fit_perOri{g}{o}{count} = OSI_fit(cell_id);
                                        tmp = OSI_global(id1{2}, 2);
                                        OSI_global_perOri{g}{o}{count} = OSI_global(cell_id);
                                    end
                                    %                                     cell_no_selective{g}(count) = length(id1{2});
                                end
                                
                            end
                            
                            
                            
                            
                            
                            % individual cell data for spatial regularity
                            if t == 2
                                %                                 for case_i = 1 : 17
                                %                                     if ~isempty(pref_diff{case_i, th})
                                %                                         individual_pdist{th, g}{case_i}{count} = [dist_p{case_i, th}(:, 1), pref_diff{case_i, th}'];
                                %                                     else
                                %                                         individual_pdist{th, g}{case_i}{count} = [];
                                %                                     end
                                %                                     if count == 1
                                %                                         total_pdist{th, g}{case_i} = individual_pdist{th, g}{case_i}{count};
                                %                                     else
                                %                                         total_pdist{th, g}{case_i} = [total_pdist{th, g}{case_i}; individual_pdist{th, g}{case_i}{count}];
                                %                                     end
                                %                                 end
                                
                                
                                for case_i = 1 : 10
                                    %                                     if ~isempty(corr_value{case_i, th})
                                    %                                         if case_i < 3
                                    %                                             individual_corr{th, g}{case_i}{count} = [dist_p{1, th}(:, 1), corr_value{case_i, th}'];
                                    %                                         elseif case_i < 5
                                    %                                             individual_corr{th, g}{case_i}{count} = [dist_p{4, th}(:, 1), corr_value{case_i, th}'];
                                    %                                         elseif case_i < 7
                                    %                                             individual_corr{th, g}{case_i}{count} = [dist_p{7, th}(:, 1), corr_value{case_i, th}'];
                                    %                                         elseif case_i == 7
                                    %                                             individual_corr{th, g}{case_i}{count} = [dist_p{13, th}(:, 1), corr_value{case_i, th}'];
                                    %                                         elseif case_i == 8
                                    %                                             individual_corr{th, g}{case_i}{count} = [dist_p{18, th}(:, 1), corr_value{case_i, th}'];
                                    %                                         elseif case_i == 9
                                    %                                             individual_corr{th, g}{case_i}{count} = [dist_p{16, th}(:, 1), corr_value{case_i, th}'];
                                    %                                         else
                                    %                                             individual_corr{th, g}{case_i}{count} = [dist_p{17, th}(:, 1), corr_value{case_i, th}'];
                                    %                                         end
                                    %
                                    %                                     else
                                    %                                         individual_corr{th, g}{case_i}{count} = [];
                                    %                                     end
                                    %                                     if count == 1
                                    %                                         total_corr{th, g}{case_i} = individual_corr{th, g}{case_i}{count};
                                    %                                     else
                                    %                                         total_corr{th, g}{case_i} = [total_corr{th, g}{case_i}; individual_corr{th, g}{case_i}{count}];
                                    %                                     end
                                end
                                
                            end
                            
                        end
                    end
                    
                    median_individual{g}(:, :, count) = median_values; % col1: DSI, col2: DSI_fit, col3: DSI_global, col4: OSI, col5: OSI_fit, col6: OSI_global
                    count = count + 1
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
    
    
    
    
    sumPlots;
    
    
    
    %% plots per FOV or per animal
%     plot_FOV_summary;
    plot_animal_summary_111220;
    
    
end









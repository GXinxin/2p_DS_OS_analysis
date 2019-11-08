% Compute direction and orientation selectivity from 2P recordings
% 18/04/16
% requires: all_vars, stimulation information and spike2 file
clear;
clc;

% fileInfo = readtext('E:\Lab\Data\2p\analysisInputs\DG_responseInfo.txt', ' ');
fileInfo = readtext('E:\Lab\Data\2p\analysisInputs\DG_responseInfo_tra2b.txt', ' ');
data_fd = 'E:\Lab\Data\2p';
for animal_id =  21:size(fileInfo, 1)%[4 5 9 11 12 13 16 25 28 33 34 36 44]
    disp(animal_id)
    disp(fileInfo{animal_id, 1})
    getResponseDG([data_fd, fileInfo{animal_id, 1}], fileInfo{animal_id, 2}, ...
        fileInfo{animal_id, 3}, fileInfo{animal_id, 4}, animal_id, fileInfo{animal_id, 5}, fileInfo{animal_id, 6}, ...
        fileInfo{animal_id, 7}, fileInfo{animal_id, 8});
end



function getResponseDG(data_path, matlab_fn, smr_fnm, data_type, animal_id, dataset_id, th, stimlengthType, spike2Type)

% load F movie, matlab file with drift grating order and the corresponding spike2 recording
cd(data_path)
load(matlab_fn)

if data_type == 2
    data_fnm = dir(fullfile('*plane*_proc*.mat'));
    load(data_fnm.name)
else
    load('all_vars.mat')
end



% initialize parameters
if spike2Type == 1
    channels = [1 3]; % th = -0.28; % settings for 333
else
    channels = [20 19]; %th = -0.24; % settings for upstair rig
end
th2 = 2; % threshold for galvo


if stimlengthType == 2
    % for 8s stim
    stimDuration = 16; % 10s grating
    stimStatic = 8; % 5s static grating
    afterStaticStim = 8; % 5s static grating after drifting grating
    totalFrames = 360; % total frames for 15s stim
else
    
    % for 5s stim
    stimDuration = 10; % 10s grating
    stimStatic = 5; % 5s static grating
    afterStaticStim = 5; % 5s static grating after drifting grating
    totalFrames = 225; % total frames for 15s stim
end


sampleFreq = 25000; % both channels recorded at 25khz



% load spike2 file
smr_fnm
fHeader = fopen(smr_fnm);
for ChannelIndex = 1:length(channels)
    [data_raw{ChannelIndex}, my_header{ChannelIndex}] = SONGetChannel(fHeader,channels(ChannelIndex),'scale');
    ChannelIndex
end



% find stimuli onsets and offsets
isStim = (smooth(data_raw{1}, 50) < th);
stimOn = find (isStim(2:end) - isStim(1:end-1) == 1) + 1;
diff = stimOn(2:end) - stimOn(1:end-1);
badId = find(diff < 180000);
stimOn = stimOn(setdiff(1:length(stimOn), badId));
length(stimOn)



% find frame onsets and offsets
isFrame = (data_raw{2} > th2);
frameStep = find (isFrame(2:end) - isFrame(1:end-1) == 1) + 1;

for s = 1 : length(stimOn)
    if ~isempty(find(frameStep > stimOn(s), 1))
        Frame_id(1, s) = find(frameStep > stimOn(s), 1); % onset of each stimulus
    end
    if ~isempty(find(frameStep > stimOn(s) + stimStatic * sampleFreq, 1))
        Frame_id(2, s) = find(frameStep > stimOn(s) + stimStatic * sampleFreq, 1); % offset of each stimulus
    end
    if ~isempty(find(frameStep > stimOn(s) + stimDuration * sampleFreq, 1))
        Frame_id(3, s) = find(frameStep > stimOn(s) + stimDuration * sampleFreq, 1); % offset of each stimulus
    end
end

Frame_id(:, 1:2)

if Frame_id(1, 1) == 1
    Frame_id = Frame_id(:, 2:end);
end
if size(Frame_id, 2) > 161 && Frame_id(1, 1) < 500
    Frame_id = Frame_id(:, 2:end);
end
if size(Frame_id, 2) > 160 && Frame_id(1, 1) < 500
    Frame_id = Frame_id(:, 2:end);
end

Frame_id(:, 1:2)

if animal_id >= 215 && animal_id <= 223
    Frame_id = Frame_id - 50;
end


% test spike2 detection
dd = Frame_id(1, 2:end) - Frame_id(1, 1:end-1);
h = figure; plot(dd);
title(['animal', num2str(animal_id), '  startFrame', num2str(Frame_id(1, 1))])
saveas(h, ['E:\Lab\Data\2p\tra2b\summary\testPlot\reorder_animal', num2str(animal_id), '.png'])
close all



% additional neuropil subtraction steps if data was processed with suite2p
if data_type == 2
    % correct intensity values for 2 glutamate animal that were
    % preprocessed initially, value was computed outside this script
    min_value = [-70, -64, -67, -68, -70, -70, -73];
    if animal_id >= 40 && animal_id <= 46
        dat.Fcell{dataset_id} = dat.Fcell{dataset_id} - min_value(animal_id-39);
        dat.FcellNeu{dataset_id} = dat.FcellNeu{dataset_id} - min_value(animal_id-39);
    end
    
    for c = 1 : size(dat.Fcell{dataset_id}, 1)
        is_cell(c) = dat.stat(c).iscell;
        neuropilCoeff(c) = dat.stat(c).neuropilCoefficient;
    end
    F_subtracted = dat.Fcell{dataset_id} - repmat(neuropilCoeff', 1, size(dat.Fcell{dataset_id}, 2)) .* dat.FcellNeu{dataset_id}; % Fcell is the raw trace of the cell, skewness (neuropil subtracted F) was computed with 0.7 coefficient in suite2p
    F_subtracted = F_subtracted';
    cell_num = length(dat.stat);
    
    add_iscell = sum(F_subtracted) > 0;
    %     iscell = iscell .* add_iscell;
    F_subtracted = F_subtracted(:, is_cell > 0);
end
no_cell = size(F_subtracted, 2);


% low pass filter
fc = 1; % Cut off frequency
fs = 15; % Sampling rate
[b,a] = butter(3, fc/(fs/2)); % Butterworth filter of order 3
padding = 50;

if stimlengthType == 2
    for cc = 1 : no_cell
        F_subtracted_filter(:, cc) = filter(b, a, [F_subtracted(padding:-1:1, cc); F_subtracted(:, cc)]); % Will be the filtered signal
    end
end
F_subtracted_filter = F_subtracted_filter(padding+1 : end, :);


% compute df from F_subtracted, corrected version since 1/15/19
for m = 1 : 2
    if m == 2
        F_subtracted = F_subtracted_filter;
    end
    hat = 500;
    se = strel('line', hat, 0);
    tmp = F_subtracted;
    for c = 1 : no_cell
        F_subtracted(:, c) = imtophat(F_subtracted(:, c)', se)'; % remove slow baseline fluctuation
        baseline = tmp(:, c) - F_subtracted(:, c);
        F_subtracted(:, c) = F_subtracted(:, c) + mean(baseline);
    end
    
    frameNum = size(F_subtracted, 1);
    f_tmp = sort(F_subtracted, 1);
    f_min = f_tmp(ceil(frameNum*0.1), :); % compute dF/F, where F is the 10th percentile of the whole F trace
    dF = F_subtracted ./ repmat(f_min, size(F_subtracted, 1), 1) - 1;
    
    
    
    % summarize into individual angles
    responseTrace = cell(1, 8);
    for angle = 1 : 8
        id = find(angles_id == angle);
        responseTrace{angle} = [];
        for j = id
            if j <= size(Frame_id, 2) && Frame_id(1, j) + totalFrames < frameNum && Frame_id(1, j) > 15
                responseTrace{angle} = cat(3, responseTrace{angle}, ...
                    dF(Frame_id(1, j) - 15 : Frame_id(1, j) + totalFrames, :));
            end
        end
    end
    max_dFF = max(dF);
    mkdir(['acq', num2str(dataset_id-1)])
    if m == 1
        save(['acq', num2str(dataset_id-1), '\responseTrace_noZscore.mat'], 'responseTrace', 'angles_id', 'Frame_id', 'max_dFF', 'dF', 'add_iscell', 'is_cell')
    else
        save(['acq', num2str(dataset_id-1), '\responseTrace_noZscore_filtered.mat'], 'responseTrace', 'angles_id', 'Frame_id', 'max_dFF', 'dF', 'add_iscell', 'is_cell')
    end
    
    
    
    %%
    % no_cell
    if m == 2
    for c = 1 : no_cell
        h = figure;
        set(h, 'position', [0 0 1200 400], 'visible', 'off')
        p = 1;
        for angle = 1 : 8
            subplot(2, 4, angle)
            for i = 1 : size(responseTrace{angle}, 3)
                if ~isempty(responseTrace{angle} * 100)
                    plot(responseTrace{angle}(:, c, i) * 100);
                end
                hold on;
            end
            if ~isempty(responseTrace{angle} * 100)
                plot(mean(responseTrace{angle}(:, c, :) * 100, 3), 'lineWidth', 1.5,...
                    'color', 'k');
                line([16 16+stimStatic*15-1], [0 0], 'color', 'b', 'lineWidth', 2.5)
                line([16+stimStatic*15 16+stimDuration*15-1], [0 0], 'color', 'r', 'lineWidth', 2.5)
                if afterStaticStim > 0
                    line([16+stimDuration*15 16+(stimDuration+afterStaticStim) * 15-1], [0 0], 'color', 'b', 'lineWidth', 2.5)
                end
                ylim([-.5, max(dF(:, c)) * 100])
                set(gca, 'Xticklabel', []); box off; %axis off
            end
            p = p + 1;
        end
        if m == 1
            saveas(h, ['acq', num2str(dataset_id-1), '\cell', num2str(c), '.png'])
        else
            saveas(h, ['acq', num2str(dataset_id-1), '\filtered_cell', num2str(c), '.png'])
        end
        
    end
    close all
    
    end
    
    
    
end


end


% modified concatAllVars to correct for shifted baseline from moco
% apply on previous data that used moco motion correction directly
% 11/13/2018


clear; clc;


% cd 'Z:\Xinxin\Ali_example_data_02152018\saline\slc-ai162-eye-injected-saline_P15_05262018'
% min_value = -192;
% cd 'Z:\Xinxin\Ali_example_data_02152018\saline\slc-ai162-eye-injected-saline_P15_06252018\Superficial'
% min_value = -107;
% cd 'Z:\Xinxin\Ali_example_data_02152018\saline\slc-ai162-eye-injection-saline_P15_2p_pup2_06252018\Superficial'
% min_value = -102; % acq2
% min_value = -129; % acq3
% cd 'Z:\Xinxin\Ali_example_data_02152018\glutamate_blocker\slc-ai162-eye-injected-ampa-blocker_P15_2p_05262018'
% min_value = -282; % acq3
% min_value = -248; % acq4
% min_value = -224; % acq5
% cd 'Z:\Xinxin\Ali_example_data_02152018\glutamate_blocker\slc-ai162-eye-injected-ampa-blocker_P15_pup2_2p_06212018\Superficial'
% min_value = -369; % acq1
% min_value = -326; % acq2
% cd 'Z:\Xinxin\Ali_example_data_02152018\glutamate_blocker\slc-ai162-eye-injection-ampa-blocker_P15_pup3_2p_06212018'
% min_value = -230;
% cd 'Z:\Xinxin\Ali_example_data_02152018\glutamate_blocker\slc-ai162-eye-injection-ampa-blocker_P14_2p_06202018\Superficial'
% min_value = -204;
% cd 'Z:\Xinxin\Ali_example_data_02152018\epibatidine\slc-ai162-eye-injected-epibatidine_P14_2p_05252018'
% min_value = -176;
% cd 'Z:\Xinxin\Ali_example_data_02152018\epibatidine\slc-ai162-eye-injected-epibatidine_P15_2p_05262018'
% min_value = -265;
% cd 'Z:\Xinxin\Ali_example_data_02152018\epibatidine\slc-ai162-eye-injected-epibatidine_P15_pup2_2p_06252018'
% min_value = -145;
% cd 'E:\Lab\Data\2p\saline\kz_ah_slc17a_ai162_saline1_P14_2p_10032018_LE\ROIs'
% min_value = -51;
% cd 'E:\Lab\Data\2p\saline\kz_ah_slc17a_ai162_saline2_P14_2p_10032018_LE\ROIs'
% min_value = -53;
% cd 'E:\Lab\Data\2p\gabazine_20x\kz_ah_slc17a_ai162_gabazine1_P15_2p_10042018_LE\ROIs'
% min_value = -66;
% cd 'E:\Lab\Data\2p\gabazine_20x\kz_ah_slc17a_ai162_gabazine1_P15_2p_10042018_RE\ROIs'
% min_value = -73;
% cd 'E:\Lab\Data\2p\gabazine_20x\kz_ah_slc17a_ai162_gabazine2_P15_2p_10042018_LE\ROIs'
% min_value = -67;
% cd 'E:\Lab\Data\2p\gabazine_20x\kz_ah_slc17a_ai162_gabazine2_P15_2p_10042018_RE\ROIs'
% min_value = -76;
% cd 'E:\Lab\Data\2p\gabazine_20x\kz_ah_slc17a_ai162_gabazine3_P15_2p_10042018_LE\ROIs'
% min_value = -88;
% cd 'E:\Lab\Data\2p\gabazine_20x\kz_ah_slc17a_ai162_gabazine3_P15_2p_10042018_RE\ROIs'
% min_value = -56;
% cd 'E:\Lab\Data\2p\gabazine_100x\kz_ah_slc17a_ai162_gabazine1_P14_2p_09132018\ROIs'
% min_value = -297;
% cd 'E:\Lab\Data\2p\gabazine_100x\kz_ah_slc17a_ai162_gabazine2_P15_2p_09142018\ROIs'
% min_value = -100;
% cd 'E:\Lab\Data\2p\gabazine_100x\kz_ah_slc17a_ai162_gabazine1_P15_2p_09252018_1\ROIs'
% min_value = -74;
% cd 'E:\Lab\Data\2p\gabazine_100x\kz_ah_slc17a_ai162_gabazine1_P15_2p_09252018_2\ROIs'
% min_value = -91;
% cd 'E:\Lab\Data\2p\gabazine_100x\kz_ah_slc17a_ai162_gabazine2_P15_2p_09252018_1\ROIs'
% min_value = -66;
% cd 'E:\Lab\Data\2p\gabazine_100x\kz_ah_slc17a_ai162_gabazine2_P15_2p_09252018_2\ROIs'
% min_value = -63;
% cd 'E:\Lab\Data\2p\glutamate\kz_ah_slc17a_ai162_glutamate1_P14_2p_09132018_rightEye\ROIs'
% min_value = -64;
% cd 'E:\Lab\Data\2p\glutamate\kz_ah_slc17a_ai162_glutamate_P14_2p_10032018_LE\ROIs'
% min_value = -96;
% cd 'E:\Lab\Data\2p\saline\kz_ah_slc17a_ai162_saline1_P14_2p_10032018_RE\ROIs'
% min_value = -78;
% cd 'E:\Lab\Data\2p\saline\kz_ah_slc17a_ai162_saline2_P14_2p_10032018_RE\ROIs'
% min_value = -64;
% fd_path = 'E:\Lab\Data\2p\gabazine_100x\kz_slc17a_ai162_gabazine2_P16_2p_03272019\kz_slc17a_ai162_gabazine2_P16_2p_03272019_ret1_2\';
% cd([fd_path, 'ROIs_retinotopy1_2'])
% fd_path = 'E:\Lab\Data\2p\gabazine_100x\kz_slc17a_ai162_gabazine1_P15_2p_03212019\kz_slc17a_ai162_gabazine1_P15_2p_03212019_DS2\';
% cd([fd_path, 'ROIs_DS_2'])
% fd_path = 'E:\Lab\Data\2p\tra2b\kz_emxtra2b_g6s_SC_P18_2p_03212019_DS_ret_1\';
% cd([fd_path, 'ROIs_ret_1'])
fd_path = 'E:\Lab\Data\2p\saline\oldData\kz_slc17a_ai162_saline1_2p_03032019\kz_slc17a_ai162_saline1_2p_03032019_DS1_3_super\';
cd([fd_path, 'ROIs_DS1_3_super'])





need_correction = 0;
r = .7; % ratio for neuropil subtraction


% filelist = dir(fullfile('ah_slc17a_ai162_noninjected_P14_2p_11082018_3_00003*'));
% filelist = dir(fullfile('ah_slc17a_ai162*'));
filelist = dir(fullfile('kz_slc17a_ai162*'));
% filelist = dir(fullfile('kz_ah_slc17a_ai162*'));
% filelist = dir(fullfile('acq2_*'));
% filelist = dir(fullfile('kz_frmd7_g6s*'));
% filelist = dir(fullfile('kz_frmd7_g6s_P16_2p_02102019_1_00007_*'));
% filelist = dir(fullfile('3rdFOV'));
% filelist = dir(fullfile('kz_emxtra2b*'));


total_rawF = [];
total_neuropil = [];
for f = 1 : length(filelist)
    
    cd(filelist(f).name)
    load('all_vars.mat')
    
    total_rawF = [total_rawF; rawF];
    total_neuropil = [total_neuropil; neuropilF]; 
    
    cd ..
end

clear F_subtracted

% correct for baseline shift from moco
if need_correction
    total_rawF = total_rawF - 2^15 - min_value;
    total_neuropil = total_neuropil - 2^15 - min_value;
end


% subtract neuropil signal (only the fluctuation not the absolute intensity values)
no_frame = size(total_rawF, 1);
for c = 1 : size(total_rawF, 2)
    tmp = sort(total_neuropil(:, c));
    % use the lowest 10% intensity value to remove the neuropil baseline to zero
    neuropil_trace = total_neuropil(:, c) - tmp(ceil(no_frame/10));   
    F_subtracted(:, c) = total_rawF(:, c) - r * neuropil_trace; % subtract neuropil fluctuation with ratio = 0.7
end


save([fd_path, 'all_vars.mat'], 'F_subtracted')


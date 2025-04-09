%%
%*------------------------------------------------------------------------*
%* Copyright (C) 2025 Aline Xavier Fidêncio.                              *
%* SPDX-License-Identifier: Apache-2.0                                    *
%*                                                                        *
%* Licensed under the Apache License, Version 2.0 (the "License");        *
%* you may not use this file except in compliance with the License.       *
%* You may obtain a copy of the License at                                *
%*                                                                        *
%* http://www.apache.org/licenses/LICENSE-2.0                             *
%*                                                                        *
%* Unless required by applicable law or agreed to in writing, software    *
%* distributed under the License is distributed on an "AS IS" BASIS,      *
%* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied*
%* See the License for the specific language governing permissions and    *
%* limitations under the License.                                         *
%*                                                                        *
%* Authors: Aline Xavier Fidêncio, Christian Klaes, Ioannis Iossifidis.   *
%*------------------------------------------------------------------------*
%--------------------------------------------------------------------------
% File:         errp_simulation.m
% Version:      0.0.4
% Last Edited:  2023-12-29
% Author:       Aline Xavier Fidêncio <aline.fidencio@gmail.com>
%               Some code from Reinmar Kobler [3]
% Description:  Script to simulate Interaction ErrP epochs using SEREEGA 
% Reference(s): 
%               [1] Krol, L. R., Pawlitzki, J., Lotte, F., Gramann, K., & Zander, T. O. (2018). 
%               SEREEGA: Simulating Event -Related EEG Activity. Journal of
%               Neuroscience Methods, 309, 13–24. doi: 10.1016/j.jneumeth.2018.08.001
%               [2] Ehinger, B. V. Simulating EEG7ERP data with SEREEGA &
%               multiple comparison corrections. (2019). 
%               https://benediktehinger.de/blog/science/simulating-data-with-sereega-
%               multiple-comparison-corrections/?cookie-state-change=1669976618130
%               [3] Kobler, Reinmar J., et al. "HEAR to remove pops and drifts: the high-variance 
%               electrode artifact removal (HEAR) algorithm." 2019 41st annual international conference 
%               of the IEEE engineering in medicine and biology society (EMBC). IEEE, 2019.
%--------------------------------------------------------------------------
% errp_simulation - Simulate Interction ErrP epochs
%   
% Inputs:
%   
%   cfg.debug - set 1 to enable debugging plots and code {default: 0}
%   t_length - total length of each simulated epoch (ms) {default: 2000}
%   prestim  - length of the prestimulus in each epoch (ms) {default: 500}  
%   srate - desired sampling rate for the simulated data (Hz) {default: 1000 Hz}
%   
%   error_rate - percentage of ErrP epochs {default: 0.2}
%   n_epochs - total number of epochs to simulate {default: 1200}
%
%   montage - EEG channel setup to use {defualt: 'S64'} 
%             use utl_get_montage('?') for help on the current available
%             options. Edit file to add new custom montage
%
% Simulation Parameters:
%   background noise: {default: brown noise from 80 random sources 25 mm
%   apart (37.5uV +/- .5 uV)
%
%   InputError condition split into fronto-centrally located P200, N250, 
%   P320, N450 (also parietal activation) as mentioned in (with dipole locations from [1]): 
%   [1] Ferrez, P. W., & Millán, J. D. R. (2008). Error-related EEG potentials 
%    generated during simulated brain–computer interaction. IEEE transactions 
%    on biomedical engineering, 55(3), 923-929. 
%   [2] Ferrez, P., & Millán, J. (2007). EEG-based brain-computer interaction: 
%   Improved accuracy by automatic single-trial error detection. Advances in 
%   neural information processing systems, 20.
%   [3] Ferrez, P. W., & Millán, J. D. R. (2005). You Are Wrong!---Automatic 
%   Detection of Interaction Errors from Brain Waves. In Proceedings of the 
%   19th international joint conference on artificial intelligence (No. CONF).
%   [4] Spüler, M., & Niethammer, C. (2015). Error-related potentials during 
%   continuous feedback: using EEG to detect errors of different type and severity. 
%   Frontiers in human neuroscience, 9, 155. (NOTE: execution ErrP)
%   [5] Omedes, J., Iturrate, I., Minguez, J., & Montesano, L. (2015). Analysis 
%   and asynchronous detection of gradually unfolding errors during monitoring 
%   tasks. Journal of neural engineering, 12(5), 056001. (NOTE: oErrP)

%   NoError condition split into P270, N450 as from findings from:
%   [3] and [5] with dipole location for N450 from [4]
%    
% Last Changed:
%   2023-12-29  axf
%               changing small details for more realistic NoError simulation
%               considering now to add latency shifts to simulate the
%               different subjects
%   2023-12-27, 28 axf
%              re-checking source locations and parameters 
%   2023-12-21 axf
%              changing some parameters and simulation details
%   2023-01-13 axf
%              added noError condition simulation
%   2023-01-06 axf
%              added range of probabilities to simulate
%              added range of snr to simulate
%   2023-01-05 axf
%              changed code structure to be more organized
%              (changed naming standard)
%              see previous version for code on finding
%              orientation (removing debug code)
%              separating noise a signal simulations
%              adjusting the random sources such that the TOTAL
%              number of components stay 64 (not 64 + n_comps I define)
%   2023-01-04 afx 
%              added freq. modulation in theta and alpha bands
close all; clear all; clc;

subjects = {'S01', 'S02', 'S03', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15',};

% subjects = {'S01',};  % to debug 

cfg = struct('debug',0); % set 1 to plot dipole orientations
cfg.show_plots = 0; % set 1 to show several plot s

for subj = 1: length(subjects)
    close all;
    clearvars -except subjects subj cfg
    
    subject = subjects{subj}; % get subject name
    fprintf('\n##\n## Simulating subject: %s\n##\n', subject); 

    save_path = strcat('../simulations/2024-05-25/', subject);             % CHANGE HERE
    if ~exist(save_path, 'dir')
        mkdir(save_path)
    end
        
    %% configuring epochs and leadfield
    t_len = 1500;                                                          % epoch length in ms
    prestim = 500;                                                         % prestimulus period in ms
    srate = 250;                                                           % sampling rate in Hz
    error_rate = 0.2;                                                      % configure error rate to be used
    n_epochs = 1200;                                                       % configure the total number of epochs to simulate
    montage = 'S64';                                                       % define EEG channel montage to use
    radius = 10;                                                           % define source radius to search 
    orientationDv = 0.2;                                                   % define the deviation to be applied to the sources orientations
    probabilities = 1;                                                     % components probabilities to simulate {default: 1, peak always happens}

    lf = lf_generate_fromnyhead('montage', montage);                       % obtaining a 64-channel lead field from ICBM-NY
    
    save([save_path '/leadfield.mat'], 'lf');                              % save the leadfield for this subject

    %% defining components
    % source locations for ERP: define center and randomly sample around it for each simulated subject
    src_ie_p200_center = lf_get_source_nearest(lf, [-1, -8, 55], region={'Brain_Right_Cingulate_Gyrus,_anterior_division', 'Brain_Left_Cingulate_Gyrus,_anterior_division'}); % BA 24 
    src_ie_p200 = lf_get_source_inradius(lf, src_ie_p200_center, radius);
    if cfg.show_plots 
        plot_source_location(src_ie_p200, lf, 'mode', '3d');
        plot_source_location(src_ie_p200_center, lf, 'mode', '3d',  'newfig', 0, 'source_colour', [0,0,0]);
        sgtitle('ErrP P200 possible sources within the given radius') 
    end
    src_ie_p200 = src_ie_p200(randi(length(src_ie_p200)));
    
    src_ie_n250_center = lf_get_source_nearest(lf, [-1, -2, 43], region={'Brain_Right_Cingulate_Gyrus,_anterior_division', 'Brain_Left_Cingulate_Gyrus,_anterior_division'}); % BA 24
    src_ie_n250 = lf_get_source_inradius(lf, src_ie_n250_center, radius);
    src_ie_n250 = src_ie_n250(randi(length(src_ie_n250)));
    
    src_ie_p320_center = lf_get_source_nearest(lf, [-1, 2, 55], region={'Brain_Left_Paracingulate_Gyrus', 'Brain_Right_Paracingulate_Gyrus'}); % BA 32
    src_ie_p320 = lf_get_source_inradius(lf, src_ie_p320_center, radius);
    src_ie_p320 = src_ie_p320(randi(length(src_ie_p320)));
    
    src_ie_n450_center = lf_get_source_nearest(lf, [-1, -14, 61], region={'Brain_Right_Juxtapositional_Lobule_Cortex_(formerly_Supplementary_Motor_Cortex)', 'Brain_Left_Juxtapositional_Lobule_Cortex_(formerly_Supplementary_Motor_Cortex)'}); % BA 6
    src_ie_n450 = lf_get_source_inradius(lf, src_ie_n450_center, radius);
    src_ie_n450 = src_ie_n450(randi(length(src_ie_n450)));
    
    src_ne_p270_center = lf_get_source_middle(lf, region = {'Brain_Right_Juxtapositional_Lobule_Cortex_(formerly_Supplementary_Motor_Cortex)', 'Brain_Left_Juxtapositional_Lobule_Cortex_(formerly_Supplementary_Motor_Cortex)'}); % middle source of BA 6 or the same as the source
    src_ne_p270 = lf_get_source_inradius(lf, src_ne_p270_center, radius);
    src_ne_p270 = src_ne_p270(randi(length(src_ne_p270)));
    
    src_ne_n450_center = lf_get_source_nearest(lf, [-6, 10, 73],  region={'Brain_Right_Juxtapositional_Lobule_Cortex_(formerly_Supplementary_Motor_Cortex)', 'Brain_Left_Juxtapositional_Lobule_Cortex_(formerly_Supplementary_Motor_Cortex)'}); % BA 6
    src_ne_n450 = lf_get_source_inradius(lf, src_ne_n450_center, radius);
    src_ne_n450 = src_ne_n450(randi(length(src_ne_n450)));
    
    src_noise = lf_get_source_spaced(lf, 80, 25);                          % select sources for background noise 
    
    for p = 1:size(probabilities,2)
        %% signal classes
        sig_ie_p200 = struct( ...
                    'type', 'erp', ...
                    'peakLatency', 200, ...
                    'peakWidth', 50, ...
                    'peakAmplitude', 20, ...                               % from Ferrez & Millán (2008), Omedes et al. (2015)
                    'probability', probabilities(p)); 
        sig_ie_p200 = utl_check_class(sig_ie_p200);
        sig_ie_p200 = utl_set_dvslope(sig_ie_p200, 'dv', .2);              % set 20% deviation

        sig_ie_n250 = struct(...
                    'type', 'erp', ...
                    'peakLatency', 250, ...
                    'peakWidth', 50, ...
                    'peakAmplitude', -40, ...                              % From Ferrez & Millán (2008), Omedes et al. (2015)
                    'probability', probabilities(p));
        sig_ie_n250 = utl_check_class(sig_ie_n250);
        sig_ie_n250 = utl_set_dvslope(sig_ie_n250, 'dv', .2);              % set 20 % deviation

        sig_ie_p320 = struct(...
                    'type', 'erp', ...
                    'peakLatency', 320, ...
                    'peakWidth', 120, ...
                    'peakAmplitude', 50, ...                               % From Ferrez & Millán (2008), Omedes et al. (2015)
                    'probability', probabilities(p));
        sig_ie_p320 = utl_check_class(sig_ie_p320);
        sig_ie_p320 = utl_set_dvslope(sig_ie_p320, 'dv', .2);              % set 20% deviation

        sig_ie_n450 = struct(...
                    'type', 'erp', ...
                    'peakLatency', 450, ...
                    'peakWidth', 150, ...
                    'peakAmplitude', -40, ...                              % From Ferrez & Millán (2008), Omedes et al. (2015)
                    'probability', probabilities(p));
        sig_ie_n450 = utl_check_class(sig_ie_n450);
        sig_ie_n450 = utl_set_dvslope(sig_ie_n450, 'dv', .2);              % set 20% deviation

        sig_ne_p270 = struct( ...
                'type', 'erp', ...
                'peakLatency', [270, 350], ...                             % overlapping peaks to be asymmetrical                                    
                'peakWidth', [200, 250], ...                               % 250 to have the slowing fadding effect
                'peakAmplitude', [20, 10], ...                             % from Ferrez & Millán (2005), Omedes et al. (2015)
                'probability', rand.* [1 1]); 
        sig_ne_p270 = utl_check_class(sig_ne_p270);
        sig_ne_p270 = utl_set_dvslope(sig_ne_p270, 'dv', .2);              % set 20 % deviation

        sig_ne_n450 = struct( ...
                    'type', 'erp', ...
                    'peakLatency', 450, ...
                    'peakWidth', 150, ...
                    'peakAmplitude', -10, ...                              % from Ferrez & Millán (2005), Omedes et al. (2015)
                    'probability', rand); 
        sig_ne_n450 = utl_check_class(sig_ne_n450);
        sig_ne_n450 = utl_set_dvslope(sig_ne_n450, 'dv', .2);              % set 20 % deviation

        sig_noise = struct(...                                             % to simulate background processes
                    'type', 'noise', ...
                    'color', 'brown',...
                    'amplitude', 37.5, ...                                 % 5uV was used by Laurens K. and Benedikt E.;  Reinmar K. uses 37.5uV
                    'amplitudeDv', .5);
        sig_noise = utl_check_class(sig_noise);

        %% creating components
        cmp_ie_p200 = utl_create_component(src_ie_p200, sig_ie_p200, lf, 'orientation', (1 + 2*orientationDv * (rand(1,3) - 0.5)) .* [0.10, 1.00, 0.63]);
        cmp_ie_n250 = utl_create_component(src_ie_n250, sig_ie_n250, lf, 'orientation', (1 + 2*orientationDv * (rand(1,3) - 0.5)) .* [0.25, 0.45, 1.00]);
        cmp_ie_p320 = utl_create_component(src_ie_p320, sig_ie_p320, lf, 'orientation', (1 + 2*orientationDv * (rand(1,3) - 0.5)) .* [0.00, -0.20, 1.00]);
        cmp_ie_n450 = utl_create_component(src_ie_n450, sig_ie_n450, lf, 'orientation', (1 + 2*orientationDv * (rand(1,3) - 0.5)) .* [0.10, 1.00, 0.28]);

        cmp_ne_p270 = utl_create_component(src_ne_p270, sig_ne_p270, lf, 'orientation', (1 + 2*orientationDv * (rand(1,3) - 0.5)) .* [0.38, 1.00, 0.91]); % fronto-centrally located
        cmp_ne_n450 = utl_create_component(src_ne_n450, sig_ne_n450, lf, 'orientation', (1 + 2*orientationDv * (rand(1,3) - 0.5)) .* [0.00, 1.00, 0.80]); % fronto-centrally located

        % apply +/-100ms latency shift to entire erp: Deviations are sampled
        % following a normal distribution with the indicated deviation being
        % the six-sigma range, capped to never exceed the indicated maximum.
        deviation = 100;
        dv = deviation / 3 * randn(); 
        if abs(dv) > deviation, dv = deviation * sign(dv); end
        
        % InputError condition
        cmp_ie_erp_signal = [cmp_ie_p200, cmp_ie_n250, cmp_ie_p320, cmp_ie_n450];
        cmp_ie_erp_signal = utl_shift_latency(cmp_ie_erp_signal, dv);
        cmp_ie_erp_signal = utl_shift_latency(cmp_ie_erp_signal, prestim);

        % NoError condition
        cmp_ne_erp_signal = [cmp_ne_p270, cmp_ne_n450];
        cmp_ne_erp_signal = utl_shift_latency(cmp_ne_erp_signal, dv);
        cmp_ne_erp_signal = utl_shift_latency(cmp_ne_erp_signal, prestim);

        % Background noise
        cmp_noise = utl_create_component(src_noise, sig_noise, lf);
        
        %% simulating data
        % InputError condition simulation
        condition = 'ie_erp';
        config = struct('n', error_rate*n_epochs, 'srate', srate, 'length', t_len, 'prestim', prestim);

        data_noise = generate_scalpdata(cmp_noise, lf, config);
        data_ie_erp_signal = generate_scalpdata(cmp_ie_erp_signal, lf, config);

        % Mixing data 
        dataset_ie_erp = data_ie_erp_signal + data_noise;
        
        % save result
        save([save_path '/' condition '_backgroundNoise.mat'], 'data_noise');
        save([save_path '/' condition '.mat'], 'dataset_ie_erp');
        
        % Add ground truth representing the component projections in the form of an ICA decomposition to a EEGLAB dataset
        eeglab; % start eeglab
        EEG = utl_create_eeglabdataset(dataset_ie_erp, config, lf, 'marker', 'InputError');
        EEG = utl_add_icaweights_toeeglabdataset(EEG, [cmp_ie_erp_signal, cmp_noise], lf);

        % save dataset, plot and save results
        plot_errp_simulation_results

        disp('INPUT ERROR CONDITION SIMULATION DONE')

        % NoError condition simulation
        condition = 'ne_erp';
        config = struct('n',(1-error_rate)*n_epochs, 'srate', srate, 'length', t_len, 'prestim', prestim);

        data_noise = generate_scalpdata(cmp_noise, lf, config);
        data_ne_erp_signal = generate_scalpdata(cmp_ne_erp_signal, lf, config);

        % Mixing data 
        dataset_ne_erp = data_ne_erp_signal + data_noise;
        
        % save result
        save([save_path '/' condition '_backgroundNoise.mat'], 'data_noise');
        save([save_path '/' condition '.mat'], 'dataset_ne_erp');
        
        % Add ground truth representing the component projections in the form of an ICA decomposition to a EEGLAB dataset
        eeglab;
        EEG = utl_create_eeglabdataset(dataset_ne_erp, config, lf, 'marker', 'NoError');
        EEG = utl_add_icaweights_toeeglabdataset(EEG, [cmp_ne_erp_signal, cmp_noise], lf);
        
        % save dataset, plot and save results
        plot_errp_simulation_results
            
        disp('NO ERROR CONDITION SIMULATION DONE')
    end
end
disp('DONE')
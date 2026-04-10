% Clear command window, clear workspace variables, and close all figures
clc, clear, close all
disp('Calculating the contrast between substantia nigra and reference region (crus cerebri)');

%% ---------------------- USER INPUT ---------------------- %%
% Define the main directory containing subject folders
input_dir = '/eru/shares/sarpal/data/SarpalLab/DARES/7T/Neuromelanin/Control';

% Use a wildcard to find all directories starting with 'sub-'
subjects  = dir(fullfile(input_dir, 'sub-*')); 

% Initialize a cell array to store the output data (12 columns)
results   = cell(1,12);

% Define the column headers for the output file
results(1,:) = {'Subject', 'NM Dir', ...
                'SN CNR Raw (mean)', 'SN CNR (mean)', ...
                'SNc CNR (mean)', 'SNr CNR (mean)', ...
                'VTA CNR (mean)', 'PBP CNR (mean)', ...
                'LC CNR (mean)', 'SN CNR (stdev)', ...
                'CC Intensity (mean)', 'CC Intensity (stdev)'};
            
% Initialize a counter to keep track of the current row being written to in 'results'
row = 1;

%% ---------------------- LOOP ---------------------- %%
% Loop through each subject found in the input directory
for s = 1:length(subjects)
    % Skip hidden directories ('.' and '..') or non-directory files
    if ~subjects(s).isdir || strcmp(subjects(s).name, '.') || strcmp(subjects(s).name, '..'), continue; end
    
    % Get the subject's name and full folder path
    subj_name = subjects(s).name;
    subj_path = fullfile(input_dir, subj_name);
    
    % Define paths where the Region of Interest (ROI) masks are stored
    mask_path = '/eru/shares/sarpal/data/SarpalLab/DARES/7T/Neuromelanin/Masks/New_Masks';
    toolbox_mask_path = '/eru/shares/sarpal/data/SarpalLab/DARES/7T/Neuromelanin/Masks/NM_toolbox_Masks';
    
    % Find all session folders for the current subject
    sessions = dir(fullfile(subj_path, 'ses-*'));
    
    % Loop through each session for the current subject
    for ss = 1:length(sessions)
        ses_name = sessions(ss).name;
        ses_path = fullfile(subj_path, ses_name);
        fprintf('\n===== Processing %s %s =====\n', subj_name, ses_name);
        
        %% -------------- Load ROI masks -------------- %%
        % Use SPM functions to read the mask files.
        % spm_vol reads the header info, and spm_read_vols reads the actual 3D image matrix.
        % 'logical()' converts the numerical mask data into boolean arrays (1s and 0s) for easy indexing.
        SN_mask  = logical(spm_read_vols(spm_vol(fullfile(mask_path, 'SN_mask_our.nii'))));
        REF_mask = logical(spm_read_vols(spm_vol(fullfile(mask_path, 'CC_mask_our.nii')))); % CC = Crus Cerebri (Reference)
        SNc_mask = logical(spm_read_vols(spm_vol(fullfile(mask_path, 'SNc_mask_our.nii'))));
        SNr_mask = logical(spm_read_vols(spm_vol(fullfile(toolbox_mask_path, 'SNr_mask_50.nii'))));
        VTA_mask = logical(spm_read_vols(spm_vol(fullfile(toolbox_mask_path, 'VTA_mask_50.nii'))));
        PBP_mask = logical(spm_read_vols(spm_vol(fullfile(toolbox_mask_path, 'PBP_mask_50.nii'))));
        LC_mask  = logical(spm_read_vols(spm_vol(fullfile(toolbox_mask_path, 'LC_mask_50.nii'))));
                
        %% -------------- Load & Smooth NM image -------------- %%
        % Construct the filename dynamically for the current subject and session
        NM_filename = sprintf('%s_%s_temp_nm.nii', subj_name, ses_name);
        NM_file     = fullfile(ses_path, NM_filename);
        
        % Check if the Neuromelanin file exists before proceeding
        if ~exist(NM_file, 'file')
            fprintf('Warning: File not found %s. Skipping...\n', NM_file);
            continue;
        end
        
        % Read the header of the Neuromelanin image
        NM_vol = spm_vol(NM_file);
        
        % Define the Full Width at Half Maximum (FWHM) for spatial smoothing in mm (x, y, z)
        FWHM   = [0.5 0.5 0.5]; 
        
        % Create the output filename for the smoothed image (prefixed with 's')
        [fpath, fname, fext] = fileparts(NM_file);
        sNM_file = fullfile(fpath, ['s' fname fext]);
        
        % Smooth the image using SPM's smoothing function
        spm_smooth(NM_vol, sNM_file, FWHM);
        
        % Read the header and 3D matrix of the newly created smoothed image
        sNM_vol = spm_vol(sNM_file);
        NM_img  = spm_read_vols(sNM_vol);

        %% -------- Reference Region Statistics (CC Mask) -------- %%
        % Extract only the voxel intensities that fall within the Reference (Crus Cerebri) mask
        REF_signal = NM_img(REF_mask);
        
        % Calculate the mode of the reference signal (requires a custom/external 'hist_mode' function)
        REF_val    = hist_mode(REF_signal); 
        
        % Calculate the standard deviation of the reference signal, ignoring NaNs
        REF_stdev  = std(REF_signal, 'omitnan');
      
        %% -------- Create CNR Map -------- %%
        % Create an output directory for the Neuromelanin CNR maps if it doesn't exist
        output_folder = fullfile(ses_path, 'NM');
        if ~exist(output_folder, 'dir'), mkdir(output_folder); end
        
        % Define the filename for the new Contrast-to-Noise Ratio (CNR) map
        CNR_fname = fullfile(output_folder, sprintf('new_CNR_%s_%s.nii', subj_name, ses_name));
        
        % Calculate Contrast-to-Noise Ratio for every voxel in the image:
        % CNR = ((Voxel Intensity - Background Intensity) / Background Intensity) * 100
        CNR_map   = (NM_img - REF_val) / REF_val * 100;
        
        % Set up the header for the new CNR map using the smoothed image's header as a template
        CNR_vol       = sNM_vol; 
        CNR_vol.fname = CNR_fname;    
        
        % Write the CNR 3D matrix to a NIfTI file on the disk
        spm_write_vol(CNR_vol, CNR_map);

        %% -------- Compute ROI CNR values -------- %%
        % Extract the CNR values for each specific Region of Interest using logical indexing
        SN_CNR  = CNR_map(SN_mask);
        SNc_CNR = CNR_map(SNc_mask);
        SNr_CNR = CNR_map(SNr_mask);
        VTA_CNR = CNR_map(VTA_mask);
        PBP_CNR = CNR_map(PBP_mask);
        LC_CNR  = CNR_map(LC_mask);
        
        % Calculate the raw mean of the Substantia Nigra (including negative values)
        SN_CNR_mean_raw = mean(SN_CNR, 'omitnan');
        
        % Calculate the means for all ROIs, strictly including ONLY positive CNR values (> 0)
        SN_CNR_mean     = mean(SN_CNR(SN_CNR > 0), 'omitnan');
        SNc_CNR_mean    = mean(SNc_CNR(SNc_CNR > 0), 'omitnan');
        SNr_CNR_mean    = mean(SNr_CNR(SNr_CNR > 0), 'omitnan');
        VTA_CNR_mean    = mean(VTA_CNR(VTA_CNR > 0), 'omitnan');
        PBP_CNR_mean    = mean(PBP_CNR(PBP_CNR > 0), 'omitnan');
        LC_CNR_mean     = mean(LC_CNR(LC_CNR > 0), 'omitnan');
        
        % Calculate the standard deviation of the Substantia Nigra CNR
        SN_CNR_stdev    = std(SN_CNR, 'omitnan');

        %% -------- Save results to cell array -------- %%
        % Increment row counter and store the computed statistics for this session
        row = row + 1;
        results(row,:) = {subj_name, ses_path, SN_CNR_mean_raw, SN_CNR_mean, ...
                          SNc_CNR_mean, SNr_CNR_mean, VTA_CNR_mean, PBP_CNR_mean, ...
                          LC_CNR_mean, SN_CNR_stdev, REF_val, REF_stdev};
    end
end

%% ---------------------- Final Export ---------------------- %%
% Define where to save the final MATLAB variable
save_path = '/eru/shares/sarpal/data/SarpalLab/DARES/7T/Neuromelanin/Control/7T_temp_RESULT_our_mask.mat';

% Save the 'results' cell array to a .mat file
save(save_path, 'results');

% Convert the cell array into a MATLAB table for easy CSV export.
% results(2:end,:) grabs the data, results(1,:) grabs the column headers.
results_table = cell2table(results(2:end,:), 'VariableNames', results(1,:));

% Write the table out to a .csv file, replacing '.mat' in the path with '.csv'
writetable(results_table, strrep(save_path, '.mat', '.csv'));

% Print completion message. Note: the filename in the text here is hardcoded.
fprintf('\n==== DONE! Results saved to PHILLIPS_temp_RESULT.csv ====\n');
% Contents of MarsBaR ROI toolbox version 0.44
%
%   mars_arm           - wrapper function for MarsBaR marmoire object
%   mars_arm_call      - services callbacks from marmoire object set functions
%   mars_armoire       - multifunction function to get/set various stores of stuff
%   mars_blob2roi      - saves ROI for cluster in xSPM structure, containing point pt
%   mars_blob_menu     - puts up ROI menu to add to SPM results interface
%   mars_blobs2rois    - creates ROIs from spm_results_ui SPM
%   mars_build_roi     - builds ROIs via the SPM GUI
%   mars_display_roi   - utility routines for display of ROIs in graphic window
%   mars_get_option    - Get option subfield as named by ``varargin``.
%   mars_image_scaling - get image scaling data for images, maybe via SPM design
%   mars_img2rois      - creates ROIs from cluster image or image containing ROIs defined by unique nos
%   mars_new_space     - make a new image space to contain image with rotations etc
%   mars_options       - options utility routines
%   mars_orthviews     - Display Orthogonal Views of a Normalized Image
%   mars_rois2img      - creates cluster or number labelled ROI image from ROIs
%   mars_struct        - multifunction function for manipulating structures
%   mars_uifile        - wrapper for matlab uiputfile/getfile; to resolve version differences
%   mars_utils         - collection of useful utility functions for MarsBaR etc
%   mars_vol_check     - FORMAT [samef, msg, chgf] = mars_vol_check(V1, V2, ...)
%   mars_vol_utils     - collection of useful utility functions for vol structs
%   marsbar            - Startup, callback and utility routine for Marsbar
%   savestruct         - saves data in structure as variables in .mat file
%
%   @mardo/add_contrasts         - method to add contrast definitions to design
%   @mardo/adjusted_data         - Return adjusted data for estimated design and contrast no
%   @mardo/betas                 - method to get estimated betas
%   @mardo/block_mean_cols       - method returns rows for means for blocks in design
%   @mardo/block_means           - method returns means for blocks in design
%   @mardo/cd_images             - method for changing path to image files in design
%   @mardo/contrasts             - method to get or set contrasts
%   @mardo/data                  - method to get or set data object 
%   @mardo/des_struct            - get/set method for des_struct field
%   @mardo/descrip               - method gets cell string description of design
%   @mardo/design_matrix         - method returns design matrix from design
%   @mardo/design_structure      - method to get or set SPM design structure
%   @mardo/display               - display method for mardo objects
%   @mardo/error_df              - method returns error df from design
%   @mardo/event_fitted          - method to compute fitted event time course
%   @mardo/event_fitted_fir      - method to compute fitted event time courses using FIR
%   @mardo/event_signal          - method to compute % signal change from fMRI events
%   @mardo/event_types           - method to get / set event types for design
%   @mardo/event_types_named     - method returns event types structures for events with same names
%   @mardo/flip_images           - flips images in design
%   @mardo/flip_option           - get/set flag for flipping images in design
%   @mardo/full_vol              - returns vol with full path, from vols, or image names
%   @mardo/get_contrast_by_name  - get named contrast(s) from design contrast structure
%   @mardo/get_contrasts         - method to get contrasts from design object
%   @mardo/get_data              - method to get data from design object
%   @mardo/get_vol_field         - method to get named field, containing or referring to vol structs
%   @mardo/has_contrasts         - method returns 1 if design has contrasts
%   @mardo/has_filter            - returns 1 if object contains filter
%   @mardo/has_images            - returns 1 if design contains images, NaN if not known
%   @mardo/has_whitener          - returns 1 if design has whitening filter
%   @mardo/image_names           - method returning image file names for design
%   @mardo/images                - method to get or set images 
%   @mardo/is_fmri               - method returns 1 if this is an fmri design
%   @mardo/is_mars_estimated     - method returns 1 if design has been estimated in MarsBaR
%   @mardo/is_marsed             - returns 1 if design has been processed with MarsBaR
%   @mardo/is_spm_estimated      - returns 1 if design has been estimated in SPM
%   @mardo/is_valid              - returns 1 if object contains valid SPM/MarsBaR design
%   @mardo/isempty               - overloaded isempty method for mardo object
%   @mardo/isfield               - method to overload isfield for mardo objects
%   @mardo/mardo                 - mardo - class constructor for MarsBaR design object
%   @mardo/mars_tag              - returns, or sets, Mars tagging structure in design
%   @mardo/marsy_data            - method to get or set marsy data
%   @mardo/masking_struct        - method to get or set SPM masking structure
%   @mardo/modality              - method returns modality of design
%   @mardo/n_effects             - get number of effects (columns) in design
%   @mardo/n_time_points         - get number of time_points in design
%   @mardo/native_vol_ver        - return string specifying native vol type
%   @mardo/paramfields           - returns struct with fields from maroi object useful for copying objects
%   @mardo/prefix_images         - method for adding or removing prefix from file names in design
%   @mardo/refresh_contrasts     - method to refresh contrasts to match design
%   @mardo/residuals             - method returns residuals from model
%   @mardo/savestruct            - saves data in def_struct as variables in .mat file
%   @mardo/set_contrasts         - method to set contrasts into design object
%   @mardo/set_data              - method to set data for design object
%   @mardo/set_vol_field         - method to set named field, containing or referring to vol structs
%   @mardo/stat_table            - gets Mars statistics and creates statistic table as cell array
%   @mardo/subsasgn              - method to overload . notation in assignments.
%   @mardo/subsref               - method to overload the . notation.
%   @mardo/summary               - method returns cell array of strings describing design
%   @mardo/swd                   - method to get/set design directory
%   @mardo/swd_writable          - returns true if swd directory can be written to 
%   @mardo/type                  - returns SPM version string corresponding to design type
%   @mardo/ui_et_edit            - method to edit invidual event types in design
%   @mardo/ui_et_edit_cb         - method to handle callbacks from ui_et_edit 
%   @mardo/ui_event_types        - ui method for selection / editing of event types
%   @mardo/ui_event_types_cb     - method to handle callbacks from ui_event_types
%   @mardo/ui_ft_design_data     - method plots FT of design and data to graphics window
%   @mardo/ui_get_contrasts      - SPM contrast UI, wrapped for MarsBaR
%   @mardo/unfiltered_efficiency - Calculate unfiltered efficiency for given SPM design and contrast
%   @mardo/verbose               - get/set method for verbose field
%
%   @mardo/private/pr_ev_diff    - method to calculate event height for % signal change
%   @mardo/private/pr_refresh_et - Refreshes data and display of event type window after edit
%   @mardo/private/pr_sort_evs   - function to sort event according to sort type
%
%   @mardo_2/add_trial_f       - method to add trial-specific F contrasts  
%   @mardo_2/apply_filter      - applies filter in design to data
%   @mardo_2/autocorr          - method to set autocorrelation types for design
%   @mardo_2/bf_dt             - method returns length of time bin for basis functions
%   @mardo_2/block_cols        - method gets design columns for block (session / subject)
%   @mardo_2/block_rows        - returns cell array of rows for each (subject/session) block
%   @mardo_2/can_mars_estimate - method returns 1 if design can be estimated in MarsBaR
%   @mardo_2/compute_contrasts - compute and return results of contrast statistics
%   @mardo_2/convert_vols      - method that converts vol structs in design and converts to format 'ver'
%   @mardo_2/design_vol        - returns vols in appropriate format for saving in design
%   @mardo_2/estimate          - estimate method - estimates GLM for SPM2 model
%   @mardo_2/event_cols        - method gets design columns for single event 
%   @mardo_2/event_onsets      - method gets onsets and durations for event/session
%   @mardo_2/event_regressor   - method gets estimated regressor for single event 
%   @mardo_2/event_specs       - method to return event specifications for all event in model
%   @mardo_2/event_x_fir       - method to return FIR design matrix columns for session
%   @mardo_2/fill              - fills missing entries from SPM FMRI design matrix 
%   @mardo_2/fwhm              - method returns FWHM, or empty if not available
%   @mardo_2/get_images        - method to get image vols from design
%   @mardo_2/has_autocorr      - returns 1 if object contains autocorrelation specification
%   @mardo_2/has_filter        - returns 1 if object contains filter
%   @mardo_2/has_images        - returns 1 if design contains images
%   @mardo_2/has_whitener      - returns 1 if design has whitening filter
%   @mardo_2/mardo_2           - class constructor for SPM2 MarsBaR design object
%   @mardo_2/mardo_99          - method to convert SPM2 design to SPM99 design
%   @mardo_2/mars_spm_graph    - Graphical display of adjusted data
%   @mardo_2/modality          - method returns modality of design
%   @mardo_2/save_spm          - method to save design as SPM format design structure
%   @mardo_2/savestruct        - saves data in def_struct into .mat file with variable name SPM
%   @mardo_2/set_images        - method to set image vols to design
%   @mardo_2/tr                - method returns TR in seconds, or empty if not available
%   @mardo_2/type              - returns SPM version string corresponding to design type
%   @mardo_2/ui_build          - method to create / fill design via GUI
%   @mardo_2/ui_get_event      - method to select an event 
%   @mardo_2/ui_get_filter     - method to get filter via GUI
%   @mardo_2/ui_report         - method for SPM2 design reporting
%   @mardo_2/ui_report_fmri    - Interactive review of fMRI design matrix
%
%   @mardo_2/private/my_design          - returns 1 if design looks like it is of SPM99 type
%   @mardo_2/private/pr_estimate        - Estimation of a General Linear Model
%   @mardo_2/private/pr_fmri_design     - MarsBaR version of spm_fMRI design - asssembles a design for fMRI studies
%   @mardo_2/private/pr_fmristat_ar     - function returns estimated AR coefficients using fmristat algorithm
%   @mardo_2/private/pr_get_filter      - gets filter using spm_fmri_spm_ui routines
%   @mardo_2/private/pr_spm_ce          - return error covariance constraints for serially correlated data
%   @mardo_2/private/pr_spm_diff        - matrix differential
%   @mardo_2/private/pr_spm_filter      - Removes low frequency confounds X0
%   @mardo_2/private/pr_spm_get_bf      - fills in basis function structure
%   @mardo_2/private/pr_spm_get_ons     - returns input [designed effects] structures
%   @mardo_2/private/pr_spm_gpdf        - Probability Density Function (PDF) of Gamma distribution
%   @mardo_2/private/pr_spm_hrf         - returns a hemodynamic response function
%   @mardo_2/private/pr_spm_orth        - recursive orthogonalization of basis functions
%   @mardo_2/private/pr_spm_q           - returns an (n x n) autocorrelation matrix for an AR(p) process
%   @mardo_2/private/pr_spm_reml        - REML estimation of covariance components from Cov{y}
%   @mardo_2/private/pr_spm_svd         - computationally efficient SVD (that can handle sparse arguments)
%   @mardo_2/private/pr_spm_ui          - MarsBaR: setting up the general linear model for independent data
%   @mardo_2/private/pr_spm_volterra    - generalized convolution of inputs (U) with basis set (bf)
%   @mardo_2/private/pr_stat_compute    - private function to compute statistics for SPM2 design
%   @mardo_2/private/pr_stat_compute_mv - private function to compute mutlivariate statistics across ROIs
%
%   @mardo_5/autocorr       - method to set autocorrelation types for design
%   @mardo_5/convert_vols   - method that converts vol structs in design and converts to format 'ver'
%   @mardo_5/estimate       - estimate method - estimates GLM for SPM2 model
%   @mardo_5/fill           - fills missing entries from SPM FMRI design matrix 
%   @mardo_5/mardo_5        - class constructor for SPM5 MarsBaR design object
%   @mardo_5/native_vol_ver - return string specifying native vol type
%   @mardo_5/type           - returns SPM version string corresponding to design type
%   @mardo_5/ui_build       - method to create / fill design via GUI
%
%   @mardo_5/private/my_design             - returns 1 if design looks like it is of SPM5 / 8 type
%   @mardo_5/private/pr_estimate           - Estimation of a General Linear Model
%   @mardo_5/private/pr_fmri_design        - MarsBaR version of spm_fMRI design - asssembles a design for fMRI studies
%   @mardo_5/private/pr_fmristat_ar        - function returns estimated AR coefficients using fmristat algorithm
%   @mardo_5/private/pr_get_filter         - gets filter using spm_fmri_spm_ui routines
%   @mardo_5/private/pr_spm_cat            - converts a cell array into a matrix
%   @mardo_5/private/pr_spm_ce             - return error covariance constraints for serially correlated data
%   @mardo_5/private/pr_spm_diff           - matrix high-order differentials
%   @mardo_5/private/pr_spm_en             - Euclidean normalization
%   @mardo_5/private/pr_spm_fileparts      - Like fileparts, but separates off a comma separated list at the end
%   @mardo_5/private/pr_spm_filter         - Removes low frequency confounds X0
%   @mardo_5/private/pr_spm_get_bf         - fills in basis function structure
%   @mardo_5/private/pr_spm_get_ons        - returns input [designed effects] structures
%   @mardo_5/private/pr_spm_gpdf           - Probability Density Function (PDF) of Gamma distribution
%   @mardo_5/private/pr_spm_hrf            - returns a hemodynamic response function
%   @mardo_5/private/pr_spm_justify        - Justify text
%   @mardo_5/private/pr_spm_logdet         - returns the log of the determinant of positive semi-definite matrix C
%   @mardo_5/private/pr_spm_non_sphericity - return error covariance constraints for basic ANOVA designs
%   @mardo_5/private/pr_spm_orth           - recursive orthogonalization of basis functions
%   @mardo_5/private/pr_spm_q              - returns an (n x n) autocorrelation matrix for an AR(p) process
%   @mardo_5/private/pr_spm_reml           - ReML estimation of covariance components from y*y'
%   @mardo_5/private/pr_spm_select         - File selector
%   @mardo_5/private/pr_spm_svd            - computationally efficient SVD (that can handle sparse arguments)
%   @mardo_5/private/pr_spm_ui             - MarsBaR: Setting up the general linear model for independent data
%   @mardo_5/private/pr_spm_unvec          - unvectorises a vectorised array
%   @mardo_5/private/pr_spm_vec            - vectorises a numeric, cell or structure array
%   @mardo_5/private/pr_spm_volterra       - generalized convolution of inputs (U) with basis set (bf)
%   @mardo_5/private/pr_stat_compute       - private function to compute statistics for SPM2 design
%   @mardo_5/private/pr_stat_compute_mv    - private function to compute mutlivariate statistics across ROIs
%
%   @mardo_99/add_trial_f       - method to add trial-specific F contrasts  
%   @mardo_99/apply_filter      - applies filter in design to data
%   @mardo_99/autocorr          - method to report lack of autocorrelation options for SPM99
%   @mardo_99/bf_dt             - method returns length of time bin for basis functions
%   @mardo_99/block_cols        - method gets design columns for block (session / subject)
%   @mardo_99/block_rows        - returns cell array of rows for each (subject/session) block
%   @mardo_99/can_mars_estimate - method returns 1 if design can be estimated in MarsBaR
%   @mardo_99/compute_contrasts - compute and return stats
%   @mardo_99/convert_vols      - method that converts vol structs in design and converts to format 'ver'
%   @mardo_99/design_vol        - returns vols in appropriate format for saving in design
%   @mardo_99/estimate          - estimate method - estimates GLM for SPM99 model
%   @mardo_99/event_cols        - method gets design columns for single event 
%   @mardo_99/event_onsets      - method gets (estimated) onsets and durations for event/session
%   @mardo_99/event_regressor   - method gets estimated regressor for single event 
%   @mardo_99/event_specs       - method to return event specifications for all event in model
%   @mardo_99/event_x_fir       - method to return FIR design matrix columns for session
%   @mardo_99/fill              - fills missing entries from SPM FMRI design matrix 
%   @mardo_99/fwhm              - method returns FWHM, or empty if not available
%   @mardo_99/get_images        - method to get image vols from design
%   @mardo_99/has_autocorr      - returns 1 if object contains autocorrelation specification
%   @mardo_99/has_filter        - returns 1 if object contains filter
%   @mardo_99/has_images        - returns 1 if design contains images
%   @mardo_99/mardo_2           - method to convert SPM2 design to SPM99 design
%   @mardo_99/mardo_99          - class constructor for SPM99 MarsBaR design object
%   @mardo_99/mars_spm_graph    - Graphical display of adjusted data
%   @mardo_99/modality          - method returns modality of design
%   @mardo_99/save_spm          - method to save design as SPM format design structure
%   @mardo_99/set_images        - method to set image vols into design
%   @mardo_99/tr                - method returns TR in seconds, or empty if not available
%   @mardo_99/type              - returns SPM version string corresponding to design type
%   @mardo_99/ui_build          - method to create / fill design via GUI
%   @mardo_99/ui_get_event      - method to select an event 
%   @mardo_99/ui_get_filter     - method to get filter via GUI
%   @mardo_99/ui_report         - mathod for SPM99 design reporting
%   @mardo_99/ui_report_fmri    - Interactive review of fMRI design matrix
%
%   @mardo_99/private/my_design          - returns 1 if design looks like it is of SPM99 type
%   @mardo_99/private/pr_estimate        - Estimation of a General Linear Model
%   @mardo_99/private/pr_fmri_design     - MarsBaR version of spm_fMRI design - asssembles a design for fMRI studies
%   @mardo_99/private/pr_get_filter      - gets filter using spm_fmri_spm_ui routines
%   @mardo_99/private/pr_spm_filter      - contruct and/or apply high and/or low pass filter
%   @mardo_99/private/pr_spm_get_bf      - creates basis functions for each trial type {i} in struct BF{i}
%   @mardo_99/private/pr_spm_get_ons     - returns onset times for events
%   @mardo_99/private/pr_spm_hrf         - returns a hemodynamic response function
%   @mardo_99/private/pr_spm_orth        - recursive orthogonalization of basis functions
%   @mardo_99/private/pr_spm_ui          - MarsBaR: setting up the general linear model for independent data
%   @mardo_99/private/pr_spm_volterra    - returns [design] matrix of explanatory variables
%   @mardo_99/private/pr_stat_compute    - calculates contrast value, stats and p values for contrasts
%   @mardo_99/private/pr_stat_compute_mv - compute multivariate statistics across ROIs
%
%   @marmoire/add_if_absent     - Adds item only if not already present
%   @marmoire/add_item          - add item to armoire
%   @marmoire/clear_item_data   - sets data for item to empty
%   @marmoire/default_item      - returns default item
%   @marmoire/do_save           - method  to save data for item
%   @marmoire/do_set            - private function to set data into item
%   @marmoire/get_item_data     - get data for item
%   @marmoire/get_item_param    - method to get item parameters
%   @marmoire/get_item_struct   - get whole item structure, including parameters
%   @marmoire/isempty_item_data - returns 1 if no data for this item
%   @marmoire/item_exists       - returns true if there is an item of this name
%   @marmoire/item_needs_save   - return 1 if item requires a save
%   @marmoire/marmoire          - marmoire - class constructor for marmoire container type
%   @marmoire/save_item_data    - save data for item to file
%   @marmoire/save_item_data_ui - save data for item to file using GUI
%   @marmoire/set_item_data     - sets data for item
%   @marmoire/set_item_data_ui  - sets data for item using GUI
%   @marmoire/set_item_param    - method to set item parameters
%   @marmoire/set_item_struct   - set whole item structure, including parameters
%   @marmoire/update_item_data  - updates data for item (sets data, flags change)
%
%   @marmoire/private/pr_is_nan     - (No help available)
%   @marmoire/private/pr_is_nix     - (No help available)
%   @marmoire/private/pr_isempty    - private function returns 1 if there is no data, or filename
%   @marmoire/private/pr_needs_save - private function returning 1 if item data needs save
%
%   @maroi/and           - overloaded add function 
%   @maroi/are_same      - returns 1 if rois are the same
%   @maroi/back2base     - back2base method - check for spacebase, transform thereto
%   @maroi/binarize      - binarize - returns / sets binarize value for object
%   @maroi/c_o_m         - c_o_m method - calculates unweighted centre of mass
%   @maroi/classdata     - classdata method - sets/gets class data
%   @maroi/descrip       - name - returns / sets name value for object
%   @maroi/display       - display - method 
%   @maroi/eq            - overloaded eq function 
%   @maroi/flip_lr       - flips ROI left / right
%   @maroi/ge            - overloaded ge function 
%   @maroi/get_marsy     - gets data in ROIs from images
%   @maroi/getdata       - getdata method - fetches time series data for ROI from images 
%   @maroi/gt            - overloaded gt (greater than) function 
%   @maroi/has_space     - has_space method - returns true if object has a native space
%   @maroi/history       - history - returns / sets history value for object
%   @maroi/label         - label - returns / sets label value for object
%   @maroi/le            - overloaded le (less than or equal to) function 
%   @maroi/loadobj       - loadobj method - fills fields needed for backwards compatibility
%   @maroi/lt            - overloaded lt (less than) function 
%   @maroi/maroi         - maroi - class constructor for umbrella ROI object
%   @maroi/maroi_matrix  - maroi_matrix method - converts roi to maroi matrix type
%   @maroi/minus         - overloaded minus function 
%   @maroi/mrdivide      - overloaded mrdivide (matrix right divide) function 
%   @maroi/mtimes        - overloaded mtimes function 
%   @maroi/native_space  - native_space method - returns native space of object
%   @maroi/ne            - overloaded ne function 
%   @maroi/not           - overloaded not function 
%   @maroi/or            - overloaded or function 
%   @maroi/paramfields   - returns struct with fields from maroi object useful for copying objects
%   @maroi/plus          - overloaded plus function 
%   @maroi/rdivide       - overloaded rdivide function 
%   @maroi/realpts       - realpts method - returns 3xN XYZ matrix in mm
%   @maroi/rle           - run length encoding method
%   @maroi/roithresh     - roithresh - returns / sets roithresh value for object
%   @maroi/save_as_image - method save_as_image - saves ROI as image
%   @maroi/save_mricro   - saves in MRIcro format
%   @maroi/saveroi       - saveroi method - checks fname, sets source field, saves object
%   @maroi/source        - source - returns / sets source value for object
%   @maroi/spm_hold      - hold - returns / sets hold value for object
%   @maroi/times         - overloaded times function 
%   @maroi/volume        - volume method - returns volume of ROI in mm
%   @maroi/xor           - overloaded xor function 
%
%   @maroi/private/my_classdata - my_classdata method - sets/gets class data
%   @maroi/private/my_loadroi   - my_loadroi function - loads ROI from file, sets source field
%   @maroi/private/my_roifname  - changes fname to appropriate ROI format
%
%   @maroi_box/centre       - centre method - sets / returns centre of ROI in mm
%   @maroi_box/flip_lr      - flips ROI left / right
%   @maroi_box/is_empty_roi - returns 1 if ROI contains no volume
%   @maroi_box/maroi_box    - maroi_box - class constructor
%   @maroi_box/volume       - volume method - returns volume of ROI in mm
%   @maroi_box/voxpts       - voxpts method - voxels within a box in given space
%
%   @maroi_image/flip_lr      - flips ROI left / right
%   @maroi_image/loadobj      - loadobj method - reloads matrix from img file
%   @maroi_image/maroi_image  - maroi_image - class constructor
%   @maroi_image/maroi_matrix - maroi_matrix method - converts maroi_image to maroi_matrix
%   @maroi_image/saveobj      - saveobj method - removes matrix information from parent to save space
%   @maroi_image/vol          - vol - returns / sets image vol for object
%
%   @maroi_image/private/my_vol_func - checks vol and func, returns processed image matrix or error
%
%   @maroi_matrix/do_write_image - method saves matrix as image and returns spm_vol
%   @maroi_matrix/domaths        - helper function to do maths on matrix object
%   @maroi_matrix/flip_lr        - flips ROI left / right
%   @maroi_matrix/is_empty_roi   - is_empty_roi - returns true if ROI contains no volume
%   @maroi_matrix/loadobj        - loadobj function - undoes run length encoding if appropriate
%   @maroi_matrix/maroi_matrix   - maroi_matrix - class constructor
%   @maroi_matrix/matrixdata     - matrixdata method - gets matrix from ROI object
%   @maroi_matrix/native_space   - native_space method - returns native space of object
%   @maroi_matrix/rebase         - rebase method - returns data from maroi_matrix in new space
%   @maroi_matrix/saveobj        - saveobj function - does run length encoding if helpful
%   @maroi_matrix/spm_mat        - spm_mat method - returns mat file defining orientation etc
%   @maroi_matrix/voxpts         - voxpts method - returns 3xN ijk matrix in voxels
%
%   @maroi_matrix/private/my_rld - function to do run length decoding 
%   @maroi_matrix/private/my_rle - method to do run length encoding on matrix
%
%   @maroi_pointlist/flip_lr         - flips ROI left / right
%   @maroi_pointlist/getvals         - returns vals for pointlist object
%   @maroi_pointlist/is_empty_roi    - is_empty_roi - returns true if ROI contains no volume
%   @maroi_pointlist/loadobj         - loadobj method - creates temporary voxel block
%   @maroi_pointlist/maroi_matrix    - maroi_matrix method - converts roi to maroi matrix type
%   @maroi_pointlist/maroi_pointlist - maroi_pointlist - class constructor
%   @maroi_pointlist/native_space    - native_space method - returns native space of object
%   @maroi_pointlist/saveobj         - saveobj method - removes temporary voxblock structure
%   @maroi_pointlist/voxpts          - voxpts method - returns 3xN ijk matrix in voxels
%
%   @maroi_pointlist/private/my_voxblock - returns voxel block and modified mat file for pointlist
%
%   @maroi_shape/c_o_m        - c_o_m method - calculates centre of mass
%   @maroi_shape/has_space    - has_space method - returns true if object has a native space
%   @maroi_shape/maroi_matrix - method to convert shape objects to maroi_matrix objects
%   @maroi_shape/maroi_shape  - maroi_shape - (virtual) shape roi class constructor
%
%   @maroi_sphere/centre       - centre method - sets / returns centre of ROI in mm
%   @maroi_sphere/flip_lr      - flips ROI left / right
%   @maroi_sphere/is_empty_roi - is_empty_roi - returns true if ROI contains no volume
%   @maroi_sphere/maroi_sphere - maroi_sphere - class constructor
%   @maroi_sphere/volume       - volume method - returns volume of ROI in mm
%   @maroi_sphere/voxpts       - voxpts method - voxels within a sphere in given space
%
%   @mars_space/display    - display method for mars_space object
%   @mars_space/eq         - overloaded eq method for mars_space objects
%   @mars_space/mars_space - mars_space - class constructor for space defining object
%   @mars_space/subsasgn   - method to over load . notation in assignments.
%   @mars_space/subsref    - method to overload the . notation.
%
%   @marsy/as_summary_only     - returns object with region data removed
%   @marsy/block_rows          - gets/sets cell array of rows for each (subject/session) block
%   @marsy/can_summarize       - returns 1 if object contains enough information to be summarized
%   @marsy/display             - display method for marsy objects
%   @marsy/eq                  - method overrides == operator
%   @marsy/is_summarized       - returns 1 if object contains calculated summary data
%   @marsy/is_summary_only     - method returns 1 if object only contains summary data
%   @marsy/is_valid            - method returns 1 if object contains valid data
%   @marsy/join                - joins marsy objects into one object
%   @marsy/marsy               - Class constructor for marsy: the MarsBaR data object
%   @marsy/n_blocks            - method returns number of subjects/sessions in data
%   @marsy/n_regions           - get number of regions
%   @marsy/n_time_points       - get number of time_points
%   @marsy/ne                  - method overrides ~= operator
%   @marsy/region              - gets / sets data for region or regions 
%   @marsy/region_data         - method gets or sets data for region(s) as cell array
%   @marsy/region_descrip      - method gets or sets descrip for region(s) as cell array
%   @marsy/region_field        - method gets or sets data for region field
%   @marsy/region_info         - method gets or sets info for region(s) as cell array
%   @marsy/region_name         - method gets or sets data for region name
%   @marsy/region_size         - method to get size of specified region data
%   @marsy/region_weights      - method gets or sets weights for region(s) as cell array
%   @marsy/resummarize         - recalculate summary data if possible
%   @marsy/savestruct          - saves data in y_struct as variables in .mat file
%   @marsy/split               - method splits regions in object into separate objects
%   @marsy/sumfunc             - method to get or set sumfunc
%   @marsy/summary             - method returns cell array of strings describing marsy object
%   @marsy/summary_block_means - return raw means over blocks in summary data
%   @marsy/summary_data        - method to get summary data, maybe set sumfunc
%   @marsy/summary_descrip     - get/set method for summary data description
%   @marsy/summary_info        - get/set method for summary data info
%   @marsy/summary_size        - method returns number of time points x number of regions
%   @marsy/ui_plot             - method plots data in various formats
%   @marsy/verbose             - get/set method for verbose field
%   @marsy/xyz                 - gets XYZ coordinates for region 
%   @marsy/y_struct            - get/set method for y_struct field
%
%   @marsy/private/pr_sum_func - creates summary stats for region data
%
%   release/make_contents - MAKECONTENTS makes Contents file, usually in current working directory.
%   release/pre_release   - Runs pre-release export, cleanup
%   release/test_rig      - runs tests on MarsBaR using specified designs
%
%   spm2/mars_blob_ui     - Displays SPM results, and ROI menu in SPM input window
%   spm2/mars_veropts     - returns SPM version specific parameters
%   spm2/spm_create_image - Wrapper for spm_create_vol, for compatibility with SPM99
%
%   spm5/mars_blob_ui     - Displays SPM results, and ROI menu in SPM input window
%   spm5/mars_veropts     - returns SPM version specific parameters
%   spm5/spm_create_image - Wrapper for spm_create_vol, for compatibility with SPM99
%   spm5/spm_get          - compatibility function to allow spm_get calls with SPM5
%
%   spm99/mars_blob_ui    - Displays SPM results, and ROI menu in SPM input window
%   spm99/mars_veropts    - returns SPM version specific parameters
%   spm99/spm_close_vol   - Close image volume - for SPM2 / SPM99 compatibility
%   spm99/spm_create_vol  - Wrapper for spm_create_image, for compatibility with SPM2
%   spm99/spm_read_hdr    - SPM2 routine to read (SPM customised) Analyze header 
%   spm99/spm_write_plane - Write a transverse plane of image data.
%   spm99/spm_write_vol   - Write an image volume to disk, setting scales and offsets as appropriate

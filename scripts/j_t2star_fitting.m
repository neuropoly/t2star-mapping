function j_t2star_fitting(opt)
% =========================================================================
% 
% T2* fitting and correct for through-slice drop out in longer echoes using
% the method of Dahnke et al. [Dahnke and Schaeffter. Magn Reson Med (2005)
% vol. 53 (5) pp. 1202-6].
% 
% INPUT
% (opt)
%   prefix				= '_top' | '_bottom'
%   path_current		= string. current path
% 
% OUTPUT
% opt
% 
%   Example
%   j_t2star_fitting('top_','.')
%
% Author: Julien Cohen-Adad <jcohen@nmr.mgh.harvard.edu>
% 2011-09-13
% 2011-09-22: seperation into three independent functions
% 2011-09-27: prefix as input (top_ of bottom_) to be used in a shell script.
% 2011-10-02: new fitting methods
% 2011-10-03: new file: smoothFreqMap. Fixed R-L orientation by copying geometric information from source data. Convolution done by Gaussian kernel (instead of box).
% 2011-11-28: fix opt. parameters that were not passed.
%
%
% =========================================================================



%% INITIALIZATION
dbstop if error % debug if error
close all


%% DEFAULT PARAMETERS

% todo
todo.computeFreqMap             	= 1;
todo.smoothFreqMap                  = 1;
todo.computeCorrectedFitting        = 1;

% File names
path_current                        = '.';
prefix                              = '';

% Frequency map
opt.fname_multiecho_magn			= [opt.prefix,'magn'];
opt.fname_multiecho_phase			= [opt.prefix,'phase'];
opt.echo_time						= (6.34:3.2:43); % in ms
opt.fname_freq						= [opt.prefix,'freq']; % freq map
opt.fname_freq_smooth				= [opt.prefix,'freq_smooth'];
% opt.fname_freq_smooth_masked		= [opt.prefix,'freq_smooth_masked'];
opt.fname_mask						= [opt.prefix,'mask']; % mask created from magnitude data. See opt.thresh_mask for threshold value.
opt.thresh_mask						= 500; % intensity under which pixels are masked. Default=500.
opt.rmse_thresh						= 0.8; % threshold above which voxels are discarded for comuting the frequency map. RMSE results from fitting the frequency slope on the phase data. Default=2.
opt.smoothType						= 'polyfit3d'; % 'gaussian' | 'box' | 'polyfit1d' | 'polyfit3d'. Default='polyfit3d'
opt.smoothKernel					= [27 27 7]; % only for 'gaussian' and 'box'
opt.smoothPolyOrder					= 3; % Default=3.
opt.smoothDownsampling				= [2 2 2]; % 3D downsample frequency map to compute gradient along Z. Default=[2 2 2].

% Gradient map
opt.fname_gradZ						= [opt.prefix,'freqGradZ']; % gradient of frequency along Z
opt.fname_gradZ_final				= [opt.prefix,'freqGradZ_final']; % gradient of frequency along Z (final version after optimization)
opt.min_length						= 6; % minimum length of values along Z, below which values are not considered. Default=4.
% opt.polyFitOrder					= 4; % max order of polynomial fit of frequency values before getting the gradient (along Z). Default=4.
opt.dz								= 1.25; % slice thickness in mm. N.B. SHOULD INCLUDE GAP!

% T2* fitting
opt.do_optimization					= 0; % 0: Just use the initial freqGradZ value - which is acceptable if nicely computed
										 % 1: Do optimization algorithm to find the final freqGradZ value, by minimizing the standard error between corrected T2* and the signal divided by the sinc term (see Danhke et al.)
fitting_method      				= 'nlls'; % 'ols': Ordinary least squares linear fit of the log of S
											  % 'gls': Generalized least squares (=weighted least squares), to respect heteroscedasticity of the residual when taking the log of S
											  % 'nlls': Levenberg-Marquardt nonlinear fitting to exponential
											  % 'num': Numerical approximation, based on the NumART2* method in [Hagberg, MRM 2002]. Tends to overestimate T2*

                                              
%% USER PARAMETERS
if ~isfield(opt,'todo'), opt.todo = todo; end
if ~isfield(opt,'prefix'), opt.prefix = prefix; end
if ~isfield(opt,'path_current'), opt.path_current = path_current; end
if ~isfield(opt,'fitting_method'), opt.fitting_method = fitting_method; end
                                              



opt.fname_t2star_uncorrected		= [opt.prefix,'t2star_uncorrected_',opt.fitting_method]; % uncorrected T2* map
opt.fname_t2star_corrected          = [opt.prefix,'t2star_corrected_',opt.fitting_method]; % corrected T2* map
opt.fname_rsquared_uncorrected		= [opt.prefix,'rsquared_uncorrected_',opt.fitting_method]; % RMSE of uncorrected fitting
opt.fname_rsquared_corrected		= [opt.prefix,'rsquared_corrected_',opt.fitting_method]; % RMSE of corrected fitting
opt.fname_iterations				= [opt.prefix,'iterations_',opt.fitting_method]; % Convergence indication (only for NLLS). If it did not converged, then it used the GLS algo.
opt.threshold_t2star_max            = 1000; % in ms. threshold T2* map (for quantization purpose when saving in NIFTI).Suggested value=1000.
% opt.convert_to_ms					= 1; % convert T2* from second to millisecond

% Misc
opt.fsloutput						= []; % If running bash: opt.fsloutput='export FSLOUTPUTTYPE=NIFTI; ', If running TCSH: opt.fsloutput='setenv FSLOUTPUTTYPE NIFTI; '. To find automatically: opt.fsloutput=[];
opt.fname_multiecho_magn_splitZ		= 'tmp.data_magn_splitZ';
opt.fname_multiecho_phase_splitZ	= 'tmp.data_phase_splitZ';
opt.fname_log						= 'log_j_t2star_fitting.txt';
opt.verbose							= 0; % for debugging purpose



% delete log file
% if exist(opt.fname_log), delete(opt.fname_log), end


%% START FUNCTION
j_disp(opt.fname_log,['\n\n\n=========================================================================================================='])
j_disp(opt.fname_log,['   Running: j_t2star_fitting'])
j_disp(opt.fname_log,['=========================================================================================================='])
j_disp(opt.fname_log,['.. Started: ',datestr(now),'\n'])

% Find which UNIX Shell is running under the Matlab Terminal
if isempty(opt.fsloutput)
	j_disp(opt.fname_log,['\nFind which UNIX Shell is running...'])
	[status result] = unix('echo $SHELL');
	if ~isempty(findstr(result,'bash'))
		j_disp(opt.fname_log,['.. Running BASH'])
		opt.fsloutput = 'export FSLOUTPUTTYPE=NIFTI; '; % if running BASH
	elseif ~isempty(findstr(result,'tcsh'))
		j_disp(opt.fname_log,['.. Running TCSH'])
		opt.fsloutput = 'setenv FSLOUTPUTTYPE NIFTI; '; % if you're running tcsh
	elseif ~isempty(findstr(result,'tsh'))
		j_disp(opt.fname_log,['.. Running TSH'])
		opt.fsloutput = 'setenv FSLOUTPUTTYPE NIFTI; '; % if you're running tcsh
    else
		j_disp(opt.fname_log,['.. Can''t find which shell is running. Using TCSH.'])
		opt.fsloutput = 'setenv FSLOUTPUTTYPE NIFTI; '; % if you're running tcsh
	end
end

% display prefix
j_disp(opt.fname_log,['\nPrefix: ',opt.prefix])

% go to the righ folder
j_disp(opt.fname_log,['\nGo to folder: ',opt.path_current])
cd(opt.path_current)

% Compute field map of frequencies from multi-echo phase data
if opt.todo.computeFreqMap
	opt = j_t2star_computeFreqMap(opt);
    if ~isstruct(opt), return; end
end

% Smooth field map of frequencies
if opt.todo.smoothFreqMap
	opt = j_t2star_smoothFreqMap(opt);
    if ~isstruct(opt), return; end
end

% Estimate corrected T2*
if opt.todo.computeCorrectedFitting
	opt = j_t2star_computeCorrectedFitting(opt);
    if ~isstruct(opt), return; end
end

% delete temporary files
delete tmp.*

%% END FUNCTION
j_disp(opt.fname_log,['\n.. Ended: ',datestr(now)])
j_disp(opt.fname_log,['==========================================================================================================\n'])

exit

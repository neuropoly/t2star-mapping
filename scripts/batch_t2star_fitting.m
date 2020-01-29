function batch_t2star_fitting(prefix,path_current)
% =========================================================================
% 
% T2* fitting and correct for through-slice drop out in longer echoes using
% the method of Dahnke et al. [Dahnke and Schaeffter. Magn Reson Med (2005)
% vol. 53 (5) pp. 1202-6].
% 
% INPUT
% prefix				= 'top_' | 'bottom_'
% path_current          = string. current path
% 
% OUTPUT
% opt
% 
%   Example
%   batch_t2star_fitting('top_','.')
%
% Author: Julien Cohen-Adad <jcohen@nmr.mgh.harvard.edu>
% 2011-09-13
% 2011-09-22: seperation into three independent functions
% 2011-09-27: prefix as input (top_ of bottom_) to be used in a shell script.
% 2011-10-02: new fitting methods
% 2011-10-03: new file: smoothFreqMap. Fixed R-L orientation by copying geometric information from source data. Convolution done by Gaussian kernel (instead of box).
% 2011-10-03: new batch file: batch_t2star_fitting to call from batch mode.
% 2011-11-28: fix opt. parameters that were not passed.
% =========================================================================


% call main function
opt.todo.computeFreqMap             	= 1;
opt.todo.smoothFreqMap                  = 1;
opt.todo.computeCorrectedFitting        = 1;
opt.prefix                              = prefix;
opt.path_current                        = path_current;
opt.fitting_method                      = 'nlls';
j_t2star_fitting(opt);

% since called in batch mode, always exit function upon completion
exit
%return;


function data_conv = j_conv3(data,kernel,opt)
% =========================================================================
% 
% Convolution of a 2D or 3D image using a 3D kernel specified as input. 
% The main difference with the Matlab function is the possibility to ignore 
% values when doing the convolution, e.g., if there are null values in the
% neighboring of a voxel.
% 
% INPUT
% data					ixjxk double
% kernel				ixjxk kernel. !!! Kernel size should be odd !!!
% (opt)					structure
%	forbidden_value		int. Value to discard when doing the convolution. (Default=0)
%	fname_log			string. File name of log file.
% 
% OUTPUT
% data_conv
% 
%   Example
%	kernel = ones(5,5,5);
%   data_conv = j_conv3(data,kernel);
%
%
% Author: Julien Cohen-Adad <jcohen@nmr.mgh.harvard.edu>
% 2011-09-23: Created
% 2011-10-03: weight by sample value instead of number of sample (in case of gaussian kernel)
% 
% =========================================================================

% PARAMETERS
verbose					= 1;

% INITIALIZATION
dbstop if error; % debug if error
if ~exist('opt'), opt = []; end
if ~isfield(opt,'forbidden_value'), opt.forbidden_value = 0; end
if ~isfield(opt,'fname_log'), opt.fname_log = 'log_j_conv3.txt'; end


% START FUNCTION
j_disp(opt.fname_log,['\n\n\n=========================================================================================================='])
j_disp(opt.fname_log,['   Running: j_conv3'])
j_disp(opt.fname_log,['=========================================================================================================='])
j_disp(opt.fname_log,['.. Started: ',datestr(now)])

% Get size of the matrix
nx = size(data,1);
ny = size(data,2);
nz = size(data,3);

% Get half size of the kernel
nkx = floor(size(kernel,1)/2);
nky = floor(size(kernel,2)/2);
nkz = floor(size(kernel,3)/2);

% pad the outer matrix with the forbidden value
data_padded = zeros(nx+2*nkx,ny+2*nky,nz+2*nkz);
data_padded(nkx+1:nx+nkx,nky+1:ny+nky,nkz+1:nz+nkz) = data;

% loop across voxels
data_conv = ones(nx,ny,nz) .* opt.forbidden_value;
nb_counts = nx*ny*nz;
icount = 1;
if verbose, j_progress(['\nPerform 3D convolution...']), end
for ix=1+nkx:nx+nkx
	for iy=1+nky:ny+nky
		for iz=1+nkz:nz+nkz
			
			% get data in kernel
			data_kernel = data_padded(ix-nkx:ix+nkx,iy-nky:iy+nky,iz-nkz:iz+nkz);

			% look which voxels to consider in the kernel neighborhood
			ind_data = find(data_kernel ~= opt.forbidden_value);
			if ~isempty(ind_data)
				% do the convolution
				data_kernel_conv = sum(data_kernel(ind_data) .* kernel(ind_data)) / sum(kernel(ind_data));
			else
				data_kernel_conv = opt.forbidden_value;
			end
			
			% fill 3D matrix
			data_conv(ix-nkx,iy-nky,iz-nkz) = data_kernel_conv;
				
			% display progress
			if verbose, j_progress(icount/nb_counts), end
			icount = icount+1;
			
		end
	end
end



% END FUNCTION
j_disp(opt.fname_log,['\n.. Ended: ',datestr(now)])
j_disp(opt.fname_log,['==========================================================================================================\n'])

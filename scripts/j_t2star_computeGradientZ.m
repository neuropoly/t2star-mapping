function opt = j_t2star_computeGradientZ(opt)
% =========================================================================
% 
% Compute map of gradient frequencies along Z.
% 
% INPUT
% opt
% 	opt.fname_freq_smooth				= 'freq_smooth.nii';
% 	opt.fname_mask						= 'mask.nii'; % mask created from magnitude data
% 	opt.fname_gradZ						= 'freqGradZ.nii'; % gradient of frequency along Z
% 
% OUTPUT
% opt
% 
%
% TODO
%
% Author: Julien Cohen-Adad <jcohen@nmr.mgh.harvard.edu>
% 2011-09-13: Created
% 2011-09-22
% 2011-09-25: create separate function from generation of freq field map
% =========================================================================



% START FUNCTION
j_disp(opt.fname_log,['\n\n\n=========================================================================================================='])
j_disp(opt.fname_log,['   Running: j_t2star_computeGradientZ'])
j_disp(opt.fname_log,['=========================================================================================================='])
j_disp(opt.fname_log,['.. Started: ',datestr(now),'\n'])


% Get dimensions of the data...
j_disp(opt.fname_log,['Get dimensions of the data...'])
cmd = ['fslsize ',opt.fname_multiecho_magn];
[status result] = unix(cmd);
if status, error(result); end
dims = j_mri_getDimensions(result);
nx = dims(1);
ny = dims(2);
nz = dims(3);
nt = dims(4);
j_disp(opt.fname_log,['.. dimension: ',num2str(nx),' x ',num2str(ny),' x ',num2str(nz),' x ',num2str(nt)])


% initialization


% Load frequency map
j_disp(opt.fname_log,['\nLoad frequency map...'])
fname_freq = [opt.fname_freq_smooth];
j_disp(opt.fname_log,['.. File name: ',fname_freq])
[img,dims,scales,bpp,endian] = read_avw(fname_freq);
freq_map_3d_smooth = squeeze(img);
clear img



% Calculate frequency gradient in the slice direction (freqGradZ)
freq_map_3d_smooth_polyfit = zeros(nx,ny,nz);
grad_z_3d = zeros(nx,ny,nz);
grad_z_3d_masked = zeros(nx,ny,nz);
icount=1;
j_progress('\nCompute frequency gradient in Z direction (slice)...')
for ix=1:nx
	for iy=1:ny
		% initialize 1D gradient values
		grad_z = zeros(1,nz);
		% get frequency along z (discard zero values)
		freq_z = squeeze(freq_map_3d_smooth(ix,iy,:));
		ind_nonzero = find(freq_z);
		if length(ind_nonzero) >= opt.min_length
			% fit to polynomial function
			p = polyfit(ind_nonzero,freq_z(ind_nonzero),opt.polyFitOrder);
			f = polyval(p,(1:nz));
			% compute frequency gradient along Z
			grad_z = gradient(f,opt.dz/1000);		
% figure, plot(freq_z(ind_nonzero),'o'), hold on, plot(f,'r'), plot(grad_z,'g')
			% fill 3D gradient matrix
			grad_z_3d(ix,iy,:) = grad_z;
		end
		% display progress
		j_progress(icount/(nx*ny))
		icount=icount+1;
	end
end


% % shift by -1*z
% j_disp(opt.fname_log,['\nShift gradient map by -1*z...'])
% grad_z_3d_shifted(:,:,2:nz) = grad_z_3d(:,:,1:nz-1);
% grad_z_3d_shifted(:,:,1) = grad_z_3d(:,:,1);
% clear grad_z_3d
% 

% Load mask
j_disp(opt.fname_log,['\nLoad mask...'])
fname = [opt.fname_mask];
j_disp(opt.fname_log,['.. File name: ',fname])
[img,dims,scales,bpp,endian] = read_avw(fname);
mask = squeeze(img);
clear img

% % Mask frequency map
% j_disp(opt.fname_log,['\nMask frequency map...'])
% grad_z_3d_masked = grad_z_3d .* mask;
% clear grad_z_3d_shifted

% Mask gradient map
j_disp(opt.fname_log,['\nMask gradient map...'])
grad_z_3d_masked = grad_z_3d .* mask;
clear grad_z_3d_shifted


% Display figure
% figure('Color','w'), imagesc(grad_z_3d_masked(:,:,5),[-10 10]), axis image, title('Gradient Frequency Map along Z (in Hz)'), colormap gray, colorbar


% Save gradient map as NIFTI file
j_disp(opt.fname_log,['Save gradient map as NIFTI file...'])
j_disp(opt.fname_log,['.. output name: ',opt.fname_gradZ])
save_avw(grad_z_3d_masked,opt.fname_gradZ,'f',scales(1:3));
j_disp(opt.fname_log,['\nCopy geometric information from ',opt.fname_multiecho_magn,'...'])
cmd = [opt.fsloutput,'fslcpgeom ',opt.fname_multiecho_magn,' ',opt.fname_gradZ,' -d'];
j_disp(opt.fname_log,['>> ',cmd]); [status result] = unix(cmd); if status, error(result); end


% END FUNCTION
j_disp(opt.fname_log,['\n.. Ended: ',datestr(now)])
j_disp(opt.fname_log,['==========================================================================================================\n'])

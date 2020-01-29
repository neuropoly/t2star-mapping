function opt = j_t2star_computeFreqMap(opt)
% =========================================================================
% 
% Compute field map of frequencies from multi echo phase data.
% 
% INPUT
% opt
%	opt.
%	opt.fname_multiecho_phase
%	opt.fname_freq
%	opt.fname_freq_smooth
%	opt.fname_freq_smooth_masked
%	opt.fname_mask
%	opt.echo_time						= (6.34:3.2:43); % in ms
% 	opt.thresh_mask						= 500; % intensity under which pixels are masked. Default=500.
% 	opt.rmse_thresh						= 0.8; % threshold above which voxels are discarded for comuting the frequency map. RMSE results from fitting the frequency slope on the phase data. Default=2.

% 
% OUTPUT
% opt
% 
%
% Author: Julien Cohen-Adad <jcohen@nmr.mgh.harvard.edu>
% 2011-09-24: Created
% 2011-10-10: Fix bug to find number of echoes.
% =========================================================================


% PARAMETERS

% INITIALIZATION
if ~exist('opt'), opt = []; end
if ~isfield(opt,'fname_log'), opt.fname_log = 'log_j_t2star_computeFreqMap.txt'; end
if ~isfield(opt,'fsloutput'), opt.fsloutput = 'export FSLOUTPUTTYPE=NIFTI; '; end


% START FUNCTION
j_disp(opt.fname_log,['\n\n\n=========================================================================================================='])
j_disp(opt.fname_log,['   Running: j_t2star_computeFreqMap'])
j_disp(opt.fname_log,['=========================================================================================================='])
j_disp(opt.fname_log,['.. Started: ',datestr(now),'\n'])


%% Check if files exist
j_disp(opt.fname_log,['\nCheck if files exist...'])
if exist([opt.fname_multiecho_magn,'.nii']) | exist([opt.fname_multiecho_magn,'.nii.gz'])
	j_disp(opt.fname_log,['.. OK [',opt.fname_multiecho_magn,']'])
else
       j_disp(opt.fname_log,['.. FILES DOES NOT EXIST!!! [',opt.fname_multiecho_magn,']'])
	opt=-1; return
end
if exist([opt.fname_multiecho_phase,'.nii']) | exist([opt.fname_multiecho_phase,'.nii.gz'])
      j_disp(opt.fname_log,['.. OK [',opt.fname_multiecho_phase,']'])
else       
	j_disp(opt.fname_log,['.. FILES DOES NOT EXIST!!! [',opt.fname_multiecho_phase,']'])
	opt=-1; return
end


%% Get dimensions of the data...
j_disp(opt.fname_log,['\nGet dimensions of the data...'])
cmd = ['fslsize ',opt.fname_multiecho_magn];
[status result] = unix(cmd);
if status, error(result); end
dims = j_mri_getDimensions(result);
nx = dims(1);
ny = dims(2);
nz = dims(3);
nt = dims(4);
j_disp(opt.fname_log,['.. dimension: ',num2str(nx),' x ',num2str(ny),' x ',num2str(nz),' x ',num2str(nt)])


%% Check echo time
j_disp(opt.fname_log,['\nCheck echo time...'])
nb_echoes = nt;
opt.echo_time = opt.echo_time(1:nt);
j_disp(opt.fname_log,['.. number of echoes: ',num2str(nb_echoes)])
j_disp(opt.fname_log,['.. echoes: ',num2str(opt.echo_time)])
% convert to seconds
echo_time_s	= opt.echo_time/1000;


%% Split volumes (because of memory issue)
j_disp(opt.fname_log,['\nSplit volume (because of memory issue)...'])
cmd = [opt.fsloutput,'fslsplit ',opt.fname_multiecho_magn,' ',opt.fname_multiecho_magn_splitZ,' -z'];
j_disp(opt.fname_log,['>> ',cmd]); [status result] = unix(cmd); if status, error(result); end
numZ = j_numbering(nz,4,0);
cmd = [opt.fsloutput,'fslsplit ',opt.fname_multiecho_phase,' ',opt.fname_multiecho_phase_splitZ,' -z'];
j_disp(opt.fname_log,['>> ',cmd]); [status result] = unix(cmd); if status, error(result); end


% Initializations
freq_map_3d = zeros(nx,ny,nz);
freq_map_3d_masked = zeros(nx,ny,nz);
grad_z_3d = zeros(nx,ny,nz);
mask_3d = zeros(nx,ny,nz);


%% Create 3D frequency map
for iz=1:nz % TODO

	j_disp(opt.fname_log,['\nSlice: ',num2str(iz),'/',num2str(nz)])
	
	% load magnitude
	j_disp(opt.fname_log,['Load magnitude data...'])
	fname_multiecho_magn_splitZ_num = [opt.fname_multiecho_magn_splitZ,numZ{iz}];
	[img,dims,scales,bpp,endian] = read_avw(fname_multiecho_magn_splitZ_num);
	data_multiecho_magn = squeeze(img);
% 	figure('Color','w'), imagesc(data_multiecho_magn(:,:,1)), axis image, title('Magnitude of 1st echo'), colormap gray, colorbar

	% Create mask from magnitude data
	j_disp(opt.fname_log,['Create mask from magnitude data...'])
	PSF = fspecial('gaussian',5,5);
	data_multiecho_magn_smooth_2d = imfilter(squeeze(data_multiecho_magn(:,:,1)),PSF,'symmetric','conv');
% 	figure('Color','w'), imagesc(data_multiecho_magn_smooth_2d), axis image, colormap gray
	ind_mask = find(data_multiecho_magn_smooth_2d>opt.thresh_mask);
	nb_pixels = length(ind_mask);
	j_disp(opt.fname_log,['.. Number of pixels: ',num2str(nb_pixels)])
	mask_2d = zeros(nx,ny);
	mask_2d(ind_mask) = 1;
	mask_3d(:,:,iz) = mask_2d;
% 	figure('Color','w'), imagesc(mask_2d), axis image, title('Mask'), colormap gray, colorbar
	
	% load phase
	j_disp(opt.fname_log,['Load phase data...'])
	fname_multiecho_phase_splitZ_num = [opt.fname_multiecho_phase_splitZ,numZ{iz}];
	[img,dims,scales,bpp,endian] = read_avw(fname_multiecho_phase_splitZ_num);
	data_multiecho_phase = squeeze(img);
	
	% convert to Radian [0,2pi), assuming max value is 4095
	j_disp(opt.fname_log,['Convert to Radian [0,2pi), assuming max value is 4095...'])
	max_phase_rad = 2*pi*(1 - 1/4096);
	data_multiecho_phase = (data_multiecho_phase./4095)*max_phase_rad;
% 	figure('Color','w'), imagesc(data_multiecho_phase(:,:,1)), axis image, title('Phase of 1st echo (in Rad)'), colormap gray, colorbar

% 	ix=346;
% 	iy=111;
	j_progress('Compute frequency map...')
	freq_map_2d = zeros(nx,ny);
	err_phase_2d = zeros(nx,ny);
	data_multiecho_magn_2d = reshape(data_multiecho_magn,nx*ny,nt);
	data_multiecho_phase_2d = reshape(data_multiecho_phase,nx*ny,nt);
	X = cat(2,echo_time_s',ones(nb_echoes,1));
	for iPix=1:nb_pixels
% j=423, i=360, iPix=find(ind_mask==(j-1)*nx+i) % N.B. i = y on fig, and j = x on fig!!!!

        data_magn_1d = data_multiecho_magn_2d(ind_mask(iPix),:);
% figure, plot(data_magn_1d), grid, title(['Magnitude, X=',num2str(j),' Y=',num2str(i)]), xlabel('Echo time (s)'); ylabel('Magnitude'); 
		data_phase_1d = data_multiecho_phase_2d(ind_mask(iPix),:);
% figure, plot(data_phase_1d), grid, title(['Phase, X=',num2str(j),' Y=',num2str(i)]), xlabel('Echo time (s)'); ylabel('Phase (rad)'); 

		% unwrap phase
		data_phase_1d_unwrapped = unwrap(data_phase_1d);

		% Linear least square fitting of y = a.X + err
		y = data_phase_1d_unwrapped';
		a = inv(X'*X)*X'*y;

        % scale phase signal
        y_scaled = y-min(y);
        y_scaled = y_scaled./max(y_scaled);
  		% Linear least square fitting of scaled phase
		a_scaled = inv(X'*X)*X'*y_scaled;
		err_phase_2d(ind_mask(iPix)) = sqrt(sum( ( y_scaled' - (a_scaled(1).*echo_time_s+a_scaled(2)) ).^2 ));
%       
%         
% 		% compute root mean squared error on phase fitting
% 		err_phase_2d(ind_mask(iPix)) = sqrt(sum( ( y' - (a(1).*echo_time_s+a(2)) ).^2 ));
		
		% Get frequency in Hertz
		freq_map_2d(ind_mask(iPix)) = a(1)/(2*pi);

 % figure('color','w'), plot(echo_time_s,data_phase_1d_unwrapped), grid, xlabel('Echo time (s)'); ylabel('Phase (rad)'); hold on, plot(echo_time_s,a(1).*echo_time_s+a(2),'r'), legend({'Raw data','Fitted Data'}), title(['Unwrapped phase, X=',num2str(j),' Y=',num2str(i),', Freq=',num2str(freq_map_2d(ind_mask(iPix))),', RMSE=',num2str(err_phase_2d(ind_mask(iPix)))])

		% display progress
		j_progress(iPix/nb_pixels)
	end % iPix
% figure('Color','w'), imagesc(freq_map_2d,[-100 100]), axis image, title('Frequency Map (in Hz)'), colormap gray, colorbar
% figure('Color','w'), imagesc(err_phase_2d), axis image, title('RMSE of frequency estimation'), colormap gray, colorbar

	% create mask from RMSE map
	mask_freq = zeros(nx,ny);
	ind_rmse = find(err_phase_2d<opt.rmse_thresh);
	mask_freq(ind_rmse) = 1;
	freq_map_2d_masked = zeros(nx,ny);
	freq_map_2d_masked(ind_rmse) = freq_map_2d(ind_rmse);
	
	% fill 3D matrix
	freq_map_3d(:,:,iz) = freq_map_2d_masked;

end
% figure('Color','w'), imagesc(freq_map_3d(:,:,5),[-100 100]), axis image, title('Frequency Map (in Hz)'), colormap gray, colorbar

% Save frequency map as NIFTI file
j_disp(opt.fname_log,['Save masked frequency map as NIFTI file...'])
j_disp(opt.fname_log,['.. output name: ',opt.fname_freq])
save_avw(freq_map_3d,opt.fname_freq,'f',scales(1:3));
j_disp(opt.fname_log,['\nCopy geometric information from ',opt.fname_multiecho_magn,'...'])
cmd = [opt.fsloutput,'fslcpgeom ',opt.fname_multiecho_magn,' ',opt.fname_freq,' -d'];
j_disp(opt.fname_log,['>> ',cmd]); [status result] = unix(cmd); if status, error(result); end

% Save mask as NIFTI file
j_disp(opt.fname_log,['Save mask as NIFTI file...'])
j_disp(opt.fname_log,['.. output name: ',opt.fname_mask])
save_avw(mask_3d,opt.fname_mask,'b',scales(1:3));
j_disp(opt.fname_log,['\nCopy geometric information from ',opt.fname_multiecho_magn,'...'])
cmd = [opt.fsloutput,'fslcpgeom ',opt.fname_multiecho_magn,' ',opt.fname_mask,' -d'];
j_disp(opt.fname_log,['>> ',cmd]); [status result] = unix(cmd); if status, error(result); end


%% END FUNCTION
j_disp(opt.fname_log,['\n.. Ended: ',datestr(now)])
j_disp(opt.fname_log,['==========================================================================================================\n'])

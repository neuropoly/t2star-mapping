function opt = j_t2star_computeCorrectedFitting(opt)
% =========================================================================
% 
% Fit T2* corrected for through-slice drop out.
% N.B. j_t2star_computeGradientZ should be run before, as it needs the
% NIFTI file gradientZ_map.nii
% 
% 
% INPUT
% opt
% 	opt.fitting_method
% 
% OUTPUT
% 
% 
%   Example
%   j_t2star_fitting
%
% TODO
%
% Author: Julien Cohen-Adad <jcohen@nmr.mgh.harvard.edu>
% 2011-09-13
% 2011-09-22
% 2011-10-02: 1) changed RMSE for R2 goodness of fit, 2) added generalized linear least squares for estimating T2*, 3) added non-linear least square fitting
% 2011-10-09: Bug relate to the fact that when data are already splitted, it takes the wrong slab. I removed this line.
% 
% =========================================================================


% PARAMETERS



% INITIALIZATION
dbstop if error % debug if error
close all
if ~exist('opt'), opt = []; end
if ~isfield(opt,'fname_log'), opt.fname_log = 'log_j_t2star_fitting.txt'; end
if ~isfield(opt,'fsloutput'), opt.fsloutput = 'export FSLOUTPUTTYPE=NIFTI; '; end
j_disp(opt.fname_log,['\n\n\n=========================================================================================================='])
j_disp(opt.fname_log,['   Running: j_t2star_computeCorrectedFitting'])
j_disp(opt.fname_log,['=========================================================================================================='])
j_disp(opt.fname_log,['.. Started: ',datestr(now),'\n'])



% T2* fitting

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


% Load gradient map
j_disp(opt.fname_log,['\nLoad gradient map...'])
[img,dims,scales,bpp,endian] = read_avw(opt.fname_gradZ);
grad_z_3d = squeeze(img);


% Load mask
j_disp(opt.fname_log,['\nLoad mask...'])
fname = [opt.fname_mask];
j_disp(opt.fname_log,['.. File name: ',fname])
[img,dims,scales,bpp,endian] = read_avw(fname);
mask = squeeze(img);
clear img


% Check echo time
j_disp(opt.fname_log,['\nCheck echo time...'])
nb_echoes = length(opt.echo_time);
j_disp(opt.fname_log,['.. number of echoes: ',num2str(nb_echoes)])
j_disp(opt.fname_log,['.. echoes (ms): ',num2str(opt.echo_time)])
echo_time = opt.echo_time;


% Check if volumes are already splitted (for debugging purpose)
%if exist([opt.fname_multiecho_magn_splitZ,'0001.nii']) | exist([opt.fname_multiecho_magn_splitZ,'0001.nii.gz'])
%	j_disp(opt.fname_log,['\nVolume already splitted...'])	
%else
	% Split volumes (because of memory issue)
	j_disp(opt.fname_log,['\nSplit volume (because of memory issue)...'])
	cmd = [opt.fsloutput,'fslsplit ',opt.fname_multiecho_magn,' ',opt.fname_multiecho_magn_splitZ,' -z'];
	j_disp(opt.fname_log,['>> ',cmd]); [status result] = unix(cmd); if status, error(result); end
%end	
numZ = j_numbering(nz,4,0);


% loop across slices
t2star_uncorr_3d = zeros(nx,ny,nz);
t2star_corr_3d = zeros(nx,ny,nz);
grad_z_final_3d = zeros(nx,ny,nz);
% adjrsquare_3d = zeros(nx,ny,nz);
rsquared_uncorr_3d = zeros(nx,ny,nz);
rsquared_corr_3d = zeros(nx,ny,nz);
iter_3d = zeros(nx,ny,nz);
numZ = j_numbering(nz,4,0);
X = cat(2,echo_time',ones(nb_echoes,1));
for iz=1:nz

	j_disp(opt.fname_log,['\nSlice: ',num2str(iz),'/',num2str(nz)])
	
	% load magnitude
	j_disp(opt.fname_log,['Load magnitude data...'])
	fname_multiecho_magn_splitZ_num = [opt.fname_multiecho_magn_splitZ,numZ{iz}];
	[img,dims,scales,bpp,endian] = read_avw(fname_multiecho_magn_splitZ_num);
	data_multiecho_magn = squeeze(img);
% 	figure('Color','w'), imagesc(data_multiecho_magn(:,:,1)), axis image, title('Magnitude of 1st echo'), colormap gray, colorbar

	% get mask indices
	ind_mask = find(mask(:,:,iz));
	nb_pixels = length(ind_mask);
	
	% initialization
	t2star_uncorr_2d = zeros(nx,ny);
	t2star_corr_2d = zeros(nx,ny);
% 	adjrsquare_2d = zeros(nx,ny);
	rsquared_uncorr_2d = zeros(nx,ny);
	rsquared_corr_2d = zeros(nx,ny);
	iter_2d = zeros(nx,ny);
	
	% loop across pixels
	if opt.do_optimization, grad_z_final_2d = zeros(nx,ny); end
	j_progress(['Fit T2* using ''',opt.fitting_method,''' method...'])
	data_multiecho_magn_2d = reshape(data_multiecho_magn,nx*ny,nt);
	grad_z_2d = reshape(grad_z_3d,nx*ny,nz);
	for iPix=1:nb_pixels
	
		% Get data magnitude in 1D
		data_magn_1d = data_multiecho_magn_2d(ind_mask(iPix),:);

		if ~isempty(find(data_magn_1d))
			% perform uncorrected T2* fit
			S = data_magn_1d;
			TE = echo_time;
			method = opt.fitting_method;
% iPix
			[T2star,S0,Sfit,Rsquared,iter] = func_t2star_fit(S,TE,method,X,nt);
			rsquared_uncorr_2d(ind_mask(iPix)) = Rsquared;
			t2star_uncorr_2d(ind_mask(iPix)) = T2star;
	% figure, plot(TE,S,'o'), hold on, plot(TE,Sfit,'r'), legend({'Raw data','Fitted Data'}), grid

			% get initial freqGradZ value from computed map
			freqGradZ_init = grad_z_2d(ind_mask(iPix),iz);

			% get final freqGradZ value
			if opt.do_optimization
				% Minimization algorithm
				[freqGradZ_final sd_err exitflag output] = fminsearch(@(delta_f) func_t2star_optimization(data_magn_1d,echo_time,delta_f,X),delta_f_init);
				freqGradZ_final_2d(ind_mask(iPix)) = freqGradZ_final;
				% disp(' ')
				% disp(['iPix          = ',num2str(iPix)])
				% disp(['Delta_f_i     = ',num2str(delta_f_init)])
				% disp(['Delta_f_f     = ',num2str(delta_f_final)])

			else
				% Just use the initial freqGradZ value - which is acceptable if nicely computed
				freqGradZ_final = freqGradZ_init;
			end

			% Correct signal by sinc function
			data_magn_1d_corr = data_magn_1d ./ abs(sinc(freqGradZ_final*echo_time/2000)); % N.B. echo time is in ms

			% perform T2* fit
			S = data_magn_1d_corr;
			TE = echo_time;
			method = opt.fitting_method;
			[T2star,S0,Sfit,Rsquared,iter] = func_t2star_fit(S,TE,method,X,nt);
			rsquared_corr_2d(ind_mask(iPix)) = Rsquared;
			t2star_corr_2d(ind_mask(iPix)) = T2star;
			iter_2d(ind_mask(iPix)) = iter;
		
		end
			
% figure, plot(TE,S,'o'), hold on, plot(TE,Sfit,'r'), legend({'Raw data','Fitted Data'}), grid


% disp(['T2* uncorrected = ',num2str(t2star_2d(ind_mask(iPix)))])
% disp(['T2* corrected = ',num2str(t2star_corr_2d(ind_mask(iPix)))])

% 		data_magn_1d_corr = data_magn_1d ./ sinc(delta_f*echo_time/2);

% 		grad_test = (-50:1:50);
% 		clear sd_err
% 		for i=1:length(grad_test)
% 			data_magn_1d_corr = data_magn_1d ./ sinc(grad_test(i)*echo_time/2);
% 			y = log(data_magn_1d_corr)';
% 			a = inv(X'*X)*X'*y;
% 			t2star_corr = -1/a(1);
% 			Sfitted = exp(a(2)-echo_time/t2star_corr);
% 			% compute error
% 			err = Sfitted - data_magn_1d_corr;
% 			sd_err(i) = std(err);
% % 			figure, plot(echo_time,data_magn_1d_corr), hold on, plot(echo_time,Sfitted,'r'), legend({'Raw data','Fitted Data'}), grid, title(['Freq gradient=',num2str(grad_test(i)),' , err=',num2str(sd_err(i))]), ylim([0 1400])
% 		end
% 		[val ind]=min(sd_err);
% 		min_grad=grad_test(ind);
%  		figure, plot(grad_test,sd_err,'*'), grid, title(['iPix=',num2str(iPix),', MinGradFreq=',num2str(min_grad)])

		% Perform linear least square fit of log(Scorr)
		% y = a.X + err
% 		y = log(data_magn_1d_corr)';
% 		a = inv(X'*X)*X'*y;
% 		t2star_corr = -1/a(1);
		
		% display progress
 		j_progress(iPix/nb_pixels)
		
	end %iPix
	
	% fill 3D T2* matrix
	t2star_uncorr_3d(:,:,iz) = t2star_uncorr_2d;
	t2star_corr_3d(:,:,iz) = t2star_corr_2d;
	if opt.do_optimization, grad_z_final_3d(:,:,iz) = grad_z_final_2d; end
% 	adjrsquare_3d(:,:,iz) = adjrsquare_2d;
	rsquared_uncorr_3d(:,:,iz) = rsquared_uncorr_2d;
	rsquared_corr_3d(:,:,iz) = rsquared_corr_2d;
	iter_3d(:,:,iz) = iter_2d;
	
end % iz

% threshold T2* map (for quantization purpose when saving in NIFTI).
t2star_uncorr_3d(find(t2star_uncorr_3d > opt.threshold_t2star_max)) = opt.threshold_t2star_max;
t2star_corr_3d = abs(t2star_corr_3d);
t2star_corr_3d(find(t2star_corr_3d > opt.threshold_t2star_max)) = opt.threshold_t2star_max;

% % convert to millisecond
% if opt.convert_to_ms
% 	j_disp(opt.fname_log,['Convert T2* to millisecond...'])
% 	t2star_uncorr_3d = t2star_uncorr_3d.*1000;
% 	t2star_corr_3d = t2star_corr_3d.*1000;
% end

% Save uncorrected t2star map as NIFTI file
j_disp(opt.fname_log,['\nSave uncorrected t2star map as NIFTI file...'])
j_disp(opt.fname_log,['.. output name: ',opt.fname_t2star_uncorrected])
save_avw(t2star_uncorr_3d,opt.fname_t2star_uncorrected,'f',scales(1:3));
j_disp(opt.fname_log,['Copy geometric information from ',opt.fname_multiecho_magn,'...'])
cmd = [opt.fsloutput,'fslcpgeom ',opt.fname_multiecho_magn,' ',opt.fname_t2star_uncorrected,' -d'];
j_disp(opt.fname_log,['>> ',cmd]); [status result] = unix(cmd); if status, error(result); end

% Save RMSE of uncorrected fitting
j_disp(opt.fname_log,['\nSave Rsquared of uncorrected fitting...'])
j_disp(opt.fname_log,['.. output name: ',opt.fname_rsquared_uncorrected])
save_avw(abs(rsquared_uncorr_3d),opt.fname_rsquared_uncorrected,'f',scales(1:3));
j_disp(opt.fname_log,['Copy geometric information from ',opt.fname_multiecho_magn,'...'])
cmd = [opt.fsloutput,'fslcpgeom ',opt.fname_multiecho_magn,' ',opt.fname_rsquared_uncorrected,' -d'];
j_disp(opt.fname_log,['>> ',cmd]); [status result] = unix(cmd); if status, error(result); end

% Save corrected t2star map as NIFTI file
j_disp(opt.fname_log,['\nSave corrected t2star map as NIFTI file...'])
j_disp(opt.fname_log,['.. output name: ',opt.fname_t2star_corrected])
save_avw(t2star_corr_3d,opt.fname_t2star_corrected,'f',scales(1:3));
j_disp(opt.fname_log,['Copy geometric information from ',opt.fname_multiecho_magn,'...'])
cmd = [opt.fsloutput,'fslcpgeom ',opt.fname_multiecho_magn,' ',opt.fname_t2star_corrected,' -d'];
j_disp(opt.fname_log,['>> ',cmd]); [status result] = unix(cmd); if status, error(result); end

% Save Rsquared of uncorrected fitting
j_disp(opt.fname_log,['\nSave Rsquared of corrected fitting...'])
j_disp(opt.fname_log,['.. output name: ',opt.fname_rsquared_corrected])
save_avw(abs(rsquared_corr_3d),opt.fname_rsquared_corrected,'f',scales(1:3));
j_disp(opt.fname_log,['Copy geometric information from ',opt.fname_multiecho_magn,'...'])
cmd = [opt.fsloutput,'fslcpgeom ',opt.fname_multiecho_magn,' ',opt.fname_rsquared_corrected,' -d'];
j_disp(opt.fname_log,['>> ',cmd]); [status result] = unix(cmd); if status, error(result); end

% Save iter flag
j_disp(opt.fname_log,['Save number of iterations (only useful for NLLS)...'])
j_disp(opt.fname_log,['.. output name: ',opt.fname_iterations])
save_avw(iter_3d,opt.fname_iterations,'f',scales(1:3));
j_disp(opt.fname_log,['Copy geometric information from ',opt.fname_multiecho_magn,'...'])
cmd = [opt.fsloutput,'fslcpgeom ',opt.fname_multiecho_magn,' ',opt.fname_iterations,' -d'];
j_disp(opt.fname_log,['>> ',cmd]); [status result] = unix(cmd); if status, error(result); end

% Save final freqGradZ
if opt.do_optimization
	j_disp(opt.fname_log,['Save final freqGradZ...'])
	j_disp(opt.fname_log,['.. output name: ',opt.fname_gradZ_final])
	save_avw(grad_z_final_3d,opt.fname_gradZ_final,'f',scales(1:3));
end


% END FUNCTION
j_disp(opt.fname_log,['\n.. Ended: ',datestr(now)])
j_disp(opt.fname_log,['==========================================================================================================\n'])





function sd_err = func_t2star_optimization(data_magn_1d,echo_time,delta_f,X)
% Optimization function

data_magn_1d_corr = data_magn_1d ./ sinc(delta_f*echo_time/2);
y = log(data_magn_1d_corr)';
a = inv(X'*X)*X'*y;
t2star_corr = -1/a(1);
% compute error
Sfitted = exp(a(2)-echo_time/t2star_corr);
err = Sfitted - data_magn_1d_corr;
sd_err = std(err);





function [T2star,S0,Sfit,Rsquared,iter] = func_t2star_fit(S,TE,method,X,nt);
% perform T2* fit
% 
% INPUT
% S
% TE
% method
% X				for ols, gls
% nt			for num

iter = 1; % number of iterations (only for NLLS). Max = 20.

switch method

case 'ols' % ordinary least squares

	% remove zeroed values (because of the log)
	nonzero = find(S);
	y = log(S(nonzero))';
	X = X(nonzero,:);
	% LLS fitting
	a = inv(X'*X)*X'*y;
	T2star = -1/a(1);
	S0 = exp(a(2));

case 'gls' % generalized least squares

	% remove zeroed values (because of the log)
	nonzero = find(S);
	y = log(S(nonzero))';
	X = X(nonzero,:);
	% LLS fitting
	V = eye(length(y)).*repmat(1./exp(y),1,length(y));
	a = inv(X'*inv(V)*X)*X'*inv(V)*y;
	T2star = -1/a(1);
	S0 = exp(a(2));

case 'nlls'  % lin_gls

	% GLS estimation to get initial parameters
	[T2star_init,S0_init] = func_t2star_fit(S,TE,'gls',X,nt);
	% Non-linear fitting
	pin = [S0_init, T2star_init]; % vector of initial parameters to be adjusted by leasqr.
	func = @(x,p) p(1)*exp(-x/p(2)); % name of function in quotes,of the form y=f(x,p)
	stol = 0.0001; % scalar tolerances on fractional improvement in ss,default stol=.0001
	niter = 20; % scalar max no. of iterations, default = 20
	wt = 1; % wt=vec(dim=1 or length(x)) of statistical weights.  These should be set to be proportional to (sqrts of var(y))^-1; (That is, the covaraince matrix of the data is assumed to be proportional to diagonal with diagonal  equal to (wt.^2)^-1.  The constant of proportionality will be estimated.),  default=1.
	dp = 0.001*ones(size(pin)); % fractional incr of p for numerical partials,default= .001*ones(size(pin))
	dfdp = 'dfdp';
	compute_covMat = 0;
	[f,p,kvg,iter] = j_stat_nlleasqr(TE',S',pin,func,stol,niter,wt,dp,dfdp,0);
	% if no convergence, use gls estimator
	if ~kvg
		S0 = S0_init;
		T2star = T2star_init;
	else
		S0 = p(1);
		T2star = p(2);
	end
	
case 'num'  % Numerical approximation based on the NumART2* method in [Hagberg, MRM 2002].

	T2star = (TE(nt)-TE(1)) * ( S(1)+S(nt)+sum(2*S(2:nt-1)) ) / (2*(nt-1)*(S(1)-S(nt)));
	S0 = S(1).*exp(TE(1)/T2star);

end

Sfit = S0 * exp(-TE/T2star);
% figure, plot(TE,S,'o'), hold on, plot(TE,Sfit,'r'), legend({'Raw data','Fitted Data'}), grid

% Compute R2 goodness of fit
SSresid = sum((S-Sfit).^2);
SStotal = (length(S)-1) * var(S);
Rsquared = 1 - SSresid/SStotal;


function  [y] = imagesc_interp1(x,sample_nrn,sample_nan)
%%  imagesc figure with fft and fftshift
[nrn,nan] = size(x);
xx = zeros(nrn*sample_nrn,nan);
for m = 1:nan
    x(:,m) = fftshift(fft(fftshift(x(:,m))));
    xx(:,m) = ifftshift(ifft(x(:,m),length(x(:,m))*sample_nrn))*sample_nrn;
end
% figure;imagesc(abs(xx));
y = zeros(nrn*sample_nrn,nan*sample_nan);
for m = 1:nrn*sample_nrn
    xx(m,:) = fftshift(fft(fftshift(xx(m,:))));
    y(m,:) = ifftshift(ifft(xx(m,:),length(xx(m,:))*sample_nan))*sample_nan;
end
% figure;imagesc(abs(y));


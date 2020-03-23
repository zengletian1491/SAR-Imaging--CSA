function a_final=Interpolate2D(x,Interp_num)
[range_num,azimuth_num] = size(x);
a_temp = zeros(range_num*(2*Interp_num+1),azimuth_num);
a_final = zeros(range_num*(2*Interp_num+1),azimuth_num*(2*Interp_num+1));
for k = 1:azimuth_num
    data_f = fftshift(fft(x(:,k)));
    data = [zeros(1,range_num*Interp_num),data_f.',zeros(1,range_num*Interp_num)];
    a_temp(:,k) = ifft(data).';
end
for k = 1:range_num*(2*Interp_num+1)
    data_f = fftshift(fft(a_temp(k,:)));
    data = [zeros(1,azimuth_num*Interp_num),data_f,zeros(1,azimuth_num*Interp_num)];
    a_final(k,:) = ifft(data);
end
clear a_temp;
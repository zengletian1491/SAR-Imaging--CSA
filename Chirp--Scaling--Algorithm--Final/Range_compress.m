function x1=Range_compress(x1,Fs,Tp,gama)
%æ‡¿Î—πÀı≥Ã–Ú

[nrn,nan]=size(x1);
fr=[-nrn/2:nrn/2-1].'/nrn*Fs;
tp_num=ceil(Tp*Fs/2)*2;
win_sref=hamming(tp_num);
win_sref=ones(tp_num,1);
win_sref=win_sref./max(win_sref);
tp1=((0:tp_num-1)'-tp_num/2)/Fs;
phi=pi*gama*(tp1.^2);
sref=exp(i*phi).*win_sref;
sref_f=zeros(nrn,1);
sref_f(floor(nrn/2)-tp_num/2+1:floor(nrn/2)+tp_num/2)=sref;
sref_f=conj((fft(sref_f))); 

parfor m=1:nan
    x1(:,m)=ifftshift(ifft((fft(x1(:,m))).*sref_f));
end


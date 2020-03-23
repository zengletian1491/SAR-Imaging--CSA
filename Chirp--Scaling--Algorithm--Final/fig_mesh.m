function fig_mesh(x,sample,dB)

[nrn,nan] = size(x);
if mod(nrn,2) ~= 0
    nrn = nrn - 1;
end
if mod(nan,2) ~= 0
    nan = nan - 1;
end
x = imagesc_interp1(x(nrn/2+[-nrn/2+1:nrn/2],...
    nan/2+[-nan/2+1:nan/2]),sample,sample);
figure,
if nargin <3
    h = mesh(abs(x));
    light('Position',[-1 1 0],'Style','infinite');
    set(h,'edgecolor',[.5 .5 .5],'facecolor',[.5 .5 .5],'edgelighting','gouraud','facelighting','gouraud');
else
    xx = 20*log10(abs(x)./max(max(abs(x))));
    x(xx<-dB) = 0;
    h = mesh(20*log10(abs(x)./max(max(abs(x)))));
    light('Position',[-1 1 0],'Style','infinite');
    set(h,'edgecolor',[.5 .5 .5],'facecolor',[.5 .5 .5],'edgelighting','gouraud','facelighting','gouraud');
    axis([1 nrn*sample 1 nan*sample -dB 0]);
end



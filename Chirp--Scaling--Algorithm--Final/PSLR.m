function [ratio] = PSLR(s)
n = length(s);
x = abs(s).^2;
[peak,position] = max(x);
for loop = position:-1:1
    if x(loop)>x(loop-1)
       x(loop) = 0;
    else
       break; 
    end
end
for loop = position+1:n
    if x(loop)>x(loop+1)
       x(loop) = 0;
    else
       break;
    end
end
side = max(x);
ratio = 10*log10(side/peak);
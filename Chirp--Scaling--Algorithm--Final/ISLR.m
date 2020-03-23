function [Ratio] = ISLR(s)
[M,N] = size(s);
if N==1
    s = s.';
else 
    s = s;
end
L = length(s);
x = abs(s).^2;
ss = x;
[peak,position] = max(x);
 for loop = position:-1:1
     if x(loop)>x(loop-1)
        x(loop) = 0;
     else
        num = loop;
        break; 
     end
end
for loop = position+1:L
    if x(loop)>x(loop+1)
       x(loop) = 0;
    else
        Num = loop;
       break; 
    end
end
y = zeros(1,L);
y(1,num:1:Num) = ss(1,num:1:Num);
Ratio = 10*log10(abs(sum(ss,2)-sum(y,2))/sum(y,2));
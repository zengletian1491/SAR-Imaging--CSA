function width = IRW(s)
L = length(s);
x = abs(s).^2;
[peak,position] = max(x);
y = x/peak;
bounds = 0.5; 
for k = position:-1:1
    if y(k)>=bounds
        half_num_1 = k;
    else
        break
    end
end
half_num_2 = half_num_1-1;
if abs(y(half_num_1)-bounds)>abs(y(half_num_2)-bounds)
    half_num = half_num_2;
else
    half_num = half_num_1;
end
for k = position:1:L
    if y(k)>=bounds
        half_num_3 = k;
    else 
        break
    end
end
half_num_4 = half_num_3+1;
if (y(half_num_3)-bounds)>(y(half_num_4)-bounds)
    Half_num = half_num_3;
else
    Half_num = half_num_4;
end
width = Half_num-half_num;

NS = zeros(1,8);
NS(1) = 1;
NS(2) = 1;
NS = vpa(NS,25);

for C = 3:8
    for i = 1:C-1
        for j = C-i:C-1
            NS(C) = NS(C) + NS(i)*NS(j);
        end
    end
end
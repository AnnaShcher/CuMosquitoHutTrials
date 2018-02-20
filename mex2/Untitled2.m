a = zeros (3, 4);
sz = size (a);
for i = 1:sz(1)
    for j = 1:sz(2)
        a(i, j) = 1;
    end
end
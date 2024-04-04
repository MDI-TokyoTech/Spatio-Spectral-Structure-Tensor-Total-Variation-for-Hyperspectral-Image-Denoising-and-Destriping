function[u] = func_PeriodicExpansionTrans(PDu)

[v, h, c, d, s1, s2] = size(PDu);
u = zeros(v,h,c,d);
for i = 1:s1
    for j = 1:s2
        u = u + circshift(PDu(:,:,:,:,i,j), [i-1, j-1, 0]);
    end
end
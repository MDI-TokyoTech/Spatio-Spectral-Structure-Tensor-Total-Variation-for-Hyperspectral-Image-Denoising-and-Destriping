% Author: Shunsuke Ono (ono@isl.titech.ac.jp)
% Last version: May. 1, 2017

function[u] = ProjL2ball(u, f, epsilon)

temp = u-f;
radius = norm(temp(:),2);
if radius > epsilon
    u = f + (epsilon/radius)*temp;
end
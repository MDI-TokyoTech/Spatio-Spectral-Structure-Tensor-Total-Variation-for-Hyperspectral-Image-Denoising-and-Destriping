% Author: Shunsuke Ono (ono@c.titech.ac.jp)

function[u] = ProjBox(u, a, b)

u(u<a) = a;
u(u>b) = b;
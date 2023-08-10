% heterogeneous initial conditions
function u0 = hettissueics(x)
c = 0.5;
d = 5;
u0 = (d-c) * rand(4,1) + c;
end
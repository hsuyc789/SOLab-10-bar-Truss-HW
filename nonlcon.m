function [c,ceq]=nonlcon(r)
[Q,stress] = TenBarTruss(r);
c(1) = (Q(3,1)^2+Q(4,1)^2)^0.5-0.02;
c(2) = max(abs(stress))-250*10^6;
ceq = [];

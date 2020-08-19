function F = M_ss( Mss )

global sigma nu betta b alppha deltta S pi tau 

c       =Mss(1);
c_H     =Mss(2);
c_F     =Mss(3);
c_Hs    =Mss(4);
r       =Mss(5);
h       =Mss(6);
w       =Mss(7);
y       =Mss(8);
k       =Mss(9);
z       =Mss(10);
f       =Mss(11);
invs    =Mss(12);

F=[
    -h^tau+w*(c*(1-b))^(-sigma);...
	-c+c_F+c_H;...
    -c_H+c_F*(1-nu)/nu;...
    -c_F+c_Hs;...
	-r+(pi/betta);...
	-y+k^alppha*h^(1-alppha);...
	-(h/k)+((1-alppha)/alppha)*(z/w)
	-w+(1-alppha)*(alppha/z)^(alppha/(1-alppha));
    -f+z+(1-deltta);
	-f+S*r/pi;...
	-invs+deltta*k;...
	-y+c_H+c_Hs+invs;];
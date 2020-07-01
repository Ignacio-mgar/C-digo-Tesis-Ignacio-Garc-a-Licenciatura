var c y tot psi_F psi_H invs k h q f pi pi_H pi_F n z w a rer r yw rw piw psiw s
    varsigma_a varsigma_c varsigma_xi varsigma_yw varsigma_pi_H varsigma_piw varsigma_f ;
 
predetermined_variables k n;

varexo varepsilon_a varepsilon_c varepsilon_xi varepsilon_yw varepsilon_pi_H varepsilon_r varepsilon_piw varepsilon_f varepsilon_rw;
      
parameters vartheta_a vartheta_c vartheta_xi vartheta_yw vartheta_pi_H vartheta_piw vartheta_f 
beta sigma b nu rho alpha delta chi omega v phi_H phi_F delta_H delta_F mu_pi mu_y eta tau k_n Xi
css c_Hss c_Fss c_Hsss rss hss wss yss kss zss fss nss invsss
sigmaw deltaw betaw  phiw tauw rhow mu_yw mu_piw bw;
    


css       =Mss(1);
c_Hss     =Mss(2);
c_Fss     =Mss(3);
c_Hsss    =Mss(4);
rss       =Mss(5);
hss       =Mss(6);
wss       =Mss(7);
yss       =Mss(8);
kss       =Mss(9);
zss       =Mss(10);
fss       =Mss(11);
nss       =Mss(12);
invsss    =Mss(13);


beta	=0.996;
sigma   =2;
b		=0.7;
nu		=0.37;
alpha	=0.33;
delta	=0.025;
chi		=2;
omega	=0.05;
v		=0.975;
phi_H	=0.5;
phi_F	=0.5;
delta_H	=0.5;
delta_F	=0.5;
Xi      =0.000742;
mu_pi	=1.5;
mu_y	=0.5;
eta		=1.5;
rho		=0.55;
tau		=1;
k_n 	=kss/nss;

vartheta_a   =0.5;
vartheta_c   =0.5;
vartheta_f   =0.5;
vartheta_pi_H=0.5;


sigmaw = 1;
deltaw = 0.5;
betaw  = 0.99;
phiw   = 0.75;
tauw   = 1;
rhow   = 0.5;
mu_yw  = 0.5;
mu_piw = 1.25;
bw     = 0.7; 

vartheta_yw   =0.5;
vartheta_piw  =0.5;
vartheta_xi   =0.5;

model(linear);  

[name= 'Euler']
(c-b*c(-1))=(c(+1)-b*c)-((1-b)/sigma)*(r-pi(+1)+varsigma_c(+1)-varsigma_c);

[name= 'Intratemporal']
tau*h=w-sigma*(c-b*c(-1))/(1-b);

[name= 'Vaciado domestico']
y*yss = css*(nu*yw+eta*nu*psi_F+nu*eta*tot*(2-nu)+(1-nu)*c)+invs*invsss;

[name= 'Funcion de produccion']
y = varsigma_a+alpha*k+(1-alpha)*h;

[name= 'Ley de evolucion del capital']
k(+1)=delta*(invs)+(1-delta)*k;

[name= 'Q de Tobin']
q=chi*delta*(invs-k);

[name= 'Prima de financiamiento externo']
f(+1) = s+r-pi(+1);

[name= 'Brecha']
s=omega*(q+k(+1)-n(+1))+varsigma_f;

[name= 'Evolución del patrimonio neto']
n(+1)/(fss*v) = (k_n)*f-(k_n-1)*(r(-1)-pi)-omega*(k_n-1)*(q(-1)+k)+n*(omega*(k_n-1)+1);

[name= 'Demanda de capital']
f = (zss/fss)*z+(1-delta)/fss*q-q(-1);

[name= 'Componentes Inflacion']
pi=pi_H+nu*(tot-tot(-1));

[name= 'CP Domestica']
(pi_H-delta_H*pi_H(-1))=beta*(pi_H(+1)-delta_H*pi_H)+((1-phi_H)*(1-phi_H*beta)/phi_H)*(psi_H+varsigma_pi_H);

[name= 'CP Importada']
(pi_F-delta_F*pi_F(-1))=beta*(pi_F(+1)-delta_F*pi_F)+((1-phi_F)*(1-phi_F*beta)/phi_F)*(psi_F);

[name= 'Regla de Taylor']
r=rho*r(-1)+(1-rho)*((mu_y*y)+(mu_pi*pi))+varepsilon_r;

[name= 'Evolucion de los activos extranjeros']
a=(a(-1)/beta)-nu*(css/yss)*(psi_F+tot)+y-(css/yss)*c+(invsss/yss)*invs;

[name= 'PMG']
z-w=h-k;

[name= 'TOT']
(tot-tot(-1))=pi_F-pi_H;

[name= 'Identidad rer']
rer=psi_F+(1-nu)*tot;

[name= 'UIP']
(r-pi(+1))=(rw-piw(+1))+(rer(+1)-rer)-Xi*a+varsigma_xi;

[name= 'MC']
psi_H=alpha*z+(1-alpha)*w+nu*tot-varsigma_a;

// Ecnomia extranjera//
[name= 'Euler foranea']
(yw-bw*yw(-1))=(yw(+1)-bw*yw)-((1-bw)/sigmaw)*(rw-piw(+1)+varsigma_yw(+1)-varsigma_yw);
[name= 'C. Phillips foranea']
(piw-deltaw*piw(-1))= betaw*(piw(+1)-deltaw*piw)+((1-phiw)*(1-phiw*betaw)/phiw)*(psiw+varsigma_piw);
[name= 'MC foraneo']
psiw = tauw*yw+sigmaw/(1-bw)*(yw-bw*yw(-1));
[name= 'Regla de Taylor foranea']
rw = rhow*rw(-1)+(1-rhow)*((mu_yw*yw)+(mu_piw*piw))+varepsilon_rw;

//-Procesos exogenos-//
varsigma_a = vartheta_a*varsigma_a(-1)+varepsilon_a;
varsigma_c = vartheta_c*varsigma_c(-1)+varepsilon_c;
varsigma_f = vartheta_f*varsigma_f(-1)+varepsilon_f;
varsigma_xi = vartheta_xi*varsigma_xi(-1)+varepsilon_xi;
varsigma_pi_H = vartheta_pi_H*varsigma_pi_H(-1)+varepsilon_pi_H;

varsigma_yw    = vartheta_yw   *varsigma_yw(-1)   +varepsilon_yw;
varsigma_piw   = vartheta_piw  *varsigma_piw(-1)  +varepsilon_piw;



end;

write_latex_dynamic_model;

resid(1);
steady;
check;  
varobs y c invs yw piw pi r rw f;


estimated_params;
omega,0.05, gamma_pdf, 0.05, 0.015;
k_n,2, gamma_pdf, 2.0,0.3;
chi,2,gamma_pdf,2,0.5;
sigma,2,gamma_pdf,2,0.2;
tau,1,gamma_pdf,1,0.1;
rho,0.55,beta_pdf,0.55,0.2;
b,0.7,beta_pdf,0.7,0.15; 
phi_H,0.5, beta_pdf, 0.5, 0.2;
phi_F,0.5, beta_pdf, 0.5, 0.2;
delta_H,0.5, beta_pdf, 0.5, 0.2;
delta_F,0.5, beta_pdf, 0.5, 0.2;
eta,1.5,gamma_pdf,1.5,0.25;
mu_y,gamma_pdf,0.5,0.4;
mu_pi,gamma_pdf,1.5,0.13;

sigmaw,1,gamma_pdf,1,0.4;
//tauw,,gamma_pdf,2,1;
rhow,0.5,beta_pdf,0.5,0.2;
bw,0.5,beta_pdf,0.5,0.15; 
deltaw,0.5, beta_pdf, 0.5, 0.2;
phiw,0.75, beta_pdf, 0.5, 0.25;
mu_yw,gamma_pdf,0.25,0.13;
mu_piw,gamma_pdf,1.5,0.13;

vartheta_a,0.5,beta_pdf,0.5,0.2;
vartheta_c,0.5,beta_pdf,0.5,0.2;
vartheta_f,0.5,beta_pdf,0.5,0.2;

vartheta_xi,0.5,beta_pdf,0.5,0.2;
vartheta_pi_H,0.5,beta_pdf,0.5,0.2;
vartheta_yw,0.5,beta_pdf,0.5,0.2;
vartheta_piw,0.5,beta_pdf,0.5,0.2;


stderr varepsilon_f, inv_gamma_pdf, 0.01, 0.1;
stderr varepsilon_r, inv_gamma_pdf, 0.01, 0.1;
stderr varepsilon_xi, inv_gamma_pdf, 0.01, 0.1;
stderr varepsilon_pi_H, inv_gamma_pdf, 0.01, 0.1;
stderr varepsilon_c, inv_gamma_pdf, 0.01, 0.1;
stderr varepsilon_a, inv_gamma_pdf, 0.01, 0.1;

stderr varepsilon_rw, inv_gamma_pdf, 0.01, 0.1;
stderr varepsilon_piw, inv_gamma_pdf, 0.01, 0.1;
stderr varepsilon_yw, inv_gamma_pdf, 0.01, 0.1;


end;

estimation(datafile=DM,xls_sheet=x2,nobs=68,first_obs=1,mh_replic=250000,mh_nblocks=2,mh_drop=0.3,mh_jscale=0.336,mode_compute=4,mode_check, graph_format = fig, bayesian_irf, irf = 20) y c r invs pi z s q k n pi_F pi_H; 
shock_decomposition y f r invs;
identification;

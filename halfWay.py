from sympy import *
from sympy.simplify.fu import TR8
half=Rational(1,2)
quarter=Rational(1,4)

g0,omegam,omegac,omegap,lmd,theta,tau,Deltam=symbols("g0,omega_m,omega_c,omega_p,lambda,theta,tau,Delta_m",cls=Symbol,positive=True)
D,mu,rho_x1,X2=symbols("D,mu,rho,x2")

# D=lmd**2*sin(theta)**2+omegap**2

c_tau_short=I*(quarter*omegac+half*Deltam-half*omegac*rho_x1+(2*omegap-mu)/(2*D)*g0**2*rho_x1**2)*tau+half*lmd*sin(theta)*tau


c_cos_2omegap_short=I*half*g0**2*lmd*sin(theta)*(D-mu*omegap)/(omegap*D**2)*rho_x1**2*cos(2*omegap*tau)

c_sin_2omegap_short=I*g0**2*(2*omegap*D+mu*(lmd**2*sin(theta)**2-omegap**2))/(4*omegap*D**2)*rho_x1**2*sin(2*omegap*tau)

c_exp_cos_short=I*g0*sqrt(2*omegam)*lmd*sin(theta)/D*rho_x1*X2*cos(omegap*tau)\
               - I*g0**2*lmd*sin(theta)/D**2*rho_x1**2*(lmd*sin(theta)*sin(2*omegap*tau)-omegap*cos(2*omegap*tau)-omegap)


c_exp_sin_short=I*g0*sqrt(2*omegam)*(mu-omegap)/D*rho_x1*X2*sin(omegap*tau)\
               - I*g0**2*(mu-omegap)/D**2*rho_x1**2*(lmd*sin(theta)-lmd*sin(theta)*cos(2*omegap*tau)-omegap*sin(2*omegap*tau))

c_exp2tau_short=I*omegam*mu/(4*lmd*sin(theta))*(X2-g0*sqrt(2/omegam)*rho_x1*(lmd*sin(theta)*sin(omegap*tau)-omegap*cos(omegap*tau))/D)**2

c_exp2tau_short_expand=I*omegam*mu/(4*lmd*sin(theta))*X2**2\
                      -I*omegam*mu*g0/(2*lmd*sin(theta))*sqrt(2/omegam)*rho_x1*X2\
                      *(lmd*sin(theta)*sin(omegap*tau)-omegap*cos(omegap*tau))/D\
                      +I*mu*g0**2/(4*lmd*sin(theta))*rho_x1**2\
                      * ((omegap**2-lmd**2*sin(theta)**2)*cos(2*omegap*tau)-2*lmd*omegap*sin(theta)*sin(2*omegap*tau))/D**2\
                      + I*mu*g0**2/(4*lmd*D*sin(theta))*rho_x1**2



c_short=c_tau_short+c_cos_2omegap_short+c_sin_2omegap_short+c_exp_cos_short+c_exp_sin_short+c_exp2tau_short_expand


c_short_expand=expand(c_short)



func1=tau

func2=cos(omegap*tau)
func3=sin(omegap*tau)

func4=cos(2*omegap*tau)
func5=sin(2*omegap*tau)


expr_expand=expand(c_short_expand)

# pprint(expr_expand)
coef1=expr_expand.coeff(func1)
coef2=expr_expand.coeff(func2)
coef3=expr_expand.coeff(func3)
coef4=expr_expand.coeff(func4)
coef5=expr_expand.coeff(func5)


G_tau=(quarter*omegac+half*Deltam-half*omegac*rho_x1+(2*omegap-mu)/(2*D)*g0**2*rho_x1**2)*tau

# G_cos_omegap=g0*sqrt(2*omegam)*(2*lmd*sin(theta)**2+omegap*mu)/(2*D*lmd*sin(theta))*rho_x1*X2*cos(omegap*tau)

G_cos_omegap=g0*sqrt(2*omegam)*(2*lmd**2*sin(theta)**2+omegap*mu)/(2*D*lmd*sin(theta))*rho_x1*X2*cos(omegap*tau)

G_sin_omegap=g0*sqrt(2*omegam)/D*(half*mu-omegap)*rho_x1*X2*sin(omegap*tau)

G_cos_2omegap=g0**2\
             * (2*D*lmd**2*sin(theta)**2+lmd**2*mu*omegap*sin(theta)**2+mu*omegap**3)/(4*D**2*lmd*omegap*sin(theta))\
             * rho_x1**2*cos(2*omegap*tau)
G_sin_2omegap=g0**2*(2*omegap*D+mu*lmd**2*sin(theta)**2-4*omegap*lmd**2*sin(theta)**2+mu*omegap**2-4*omegap**3)/(4*omegap*D**2)\
             * rho_x1**2*sin(2*omegap*tau)

G_1=g0**2\
    *(8*lmd**2*omegap*sin(theta)**2-4*lmd**2*mu*sin(theta)**2+D*mu)/(4*lmd*sin(theta)*D**2)\
    * rho_x1**2\
    + omegam*mu/(4*lmd*sin(theta))*X2**2

tmp=expr_expand-I*(G_tau+G_cos_omegap+G_sin_omegap+G_cos_2omegap+G_sin_2omegap+G_1)



pprint(tmp.simplify())
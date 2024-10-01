from sympy import *
import sympy as sp


g0,omegam,omegac,omegap,lmd,theta,tau,Deltam=symbols("g0,omega_m,omega_c,omega_p,lambda,theta,tau,Delta_m",cls=Symbol,positive=True)


x1,xi,s2=symbols("x1,xi,s2",cls=Symbol,real=True)

half=Rational(1,2)
quarter=Rational(1,4)

#####################################################
# in this part, x2 is a function of s2
x2=g0*sqrt(2/omegam)*(omegac*x1**2-half)*(lmd*sin(theta)*sin(omegap*tau)-omegap*cos(omegap*tau))/(lmd**2*sin(theta)**2+omegap**2)\
    +(s2+g0*sqrt(2/omegam)*(omegac*x1**2-half)*omegap/(lmd**2*sin(theta)**2+omegap**2))*exp(-lmd*sin(theta)*tau)





#
#
# x2Complex=g0*sqrt(2/omegam)*(omegac*x1**2-half)*(lmd*sin(theta)*sin(omegap*tau)-omegap*cos(omegap*tau))/(lmd**2*sin(theta)**2+omegap**2)\
#     + (s2+g0*sqrt(2/omegam)*(omegac*x1**2-half)*omegap/(lmd**2*sin(theta)**2+omegap**2))*exp(-lmd*sin(theta)*tau)
#





c0=I*half*omegac+I*half*Deltam-I*half*omegac**2*x1**2+half*lmd*sin(theta)\
    +(I*half*g0*sqrt(2*omegam)-I*g0*omegac*sqrt(2*omegam)*x1**2)*x2*cos(omegap*tau)\
    -I*(half*lmd*omegam*cos(theta)+half*Deltam*omegam)*x2**2


# c0=I*half*omegac+I*half*Deltam-I*half*omegac**2*x1**2+half*lmd*sin(theta)\
#     +(I*half*g0*sqrt(2*omegam)-I*g0*omegac*sqrt(2*omegam)*x1**2)*x2Complex*cos(omegap*tau)\
#     -I*(half*lmd*omegam*cos(theta)+half*Deltam*omegam)*x2Complex**2


c_tau=I*(half*omegac+half*Deltam-half*omegac**2*x1**2+g0**2*omegap/(lmd**2*sin(theta)**2+omegap**2)*(omegac*x1**2-half)**2-half*g0**2*(lmd*cos(theta)+Deltam)/(lmd**2*sin(theta)**2+omegap**2)*(omegac*x1**2-half)**2)*tau\
    +half*lmd*sin(theta)*tau


c_sin_2omegap=I*half*g0**2*1/(lmd**2*sin(theta)**2+omegap**2)*(omegac*x1**2-half)**2*sin(2*omegap*tau)\
             + I*g0**2/(4*omegap)*(lmd*cos(theta)+Deltam)*(lmd**2*sin(theta)**2-omegap**2)/(lmd**2*sin(theta)**2+omegap**2)**2\
             * (omegac*x1**2-half)**2*sin(2*omegap*tau)

c_cos_2omegap=I*g0**2/(2*omegap)*lmd*sin(theta)/(lmd**2*sin(theta)**2+omegap**2)*(omegac*x1**2-half)**2*cos(2*omegap*tau)\
    -I*g0**2*lmd*sin(theta)/2*(lmd*cos(theta)+Deltam)/(lmd**2*sin(theta)**2+omegap**2)**2*(omegac*x1**2-half)**2*cos(2*omegap*tau)



c_exp_sin=-I*g0*omegap*sqrt(2*omegam)*(omegac*x1**2-half)\
        *(s2+g0*sqrt(2/omegam)*(omegac*x1**2-half)*omegap/(lmd**2*sin(theta)**2+omegap**2))\
        *1/(lmd**2*sin(theta)**2+omegap**2)*exp(-lmd*sin(theta)*tau)*sin(omegap*tau)\
        + I*g0*omegam*(lmd*cos(theta)+Deltam)*sqrt(2/omegam)*(omegac*x1**2-half)\
        *(s2+g0*sqrt(2/omegam)*(omegac*x1**2-half)*omegap/(lmd**2*sin(theta)**2+omegap**2))\
        *1/(lmd**2*sin(theta)**2+omegap**2)*exp(-lmd*sin(theta)*tau)*sin(omegap*tau)



c_exp_cos=I*g0*sqrt(2*omegam)*(omegac*x1**2-half)\
        * (s2+g0*sqrt(2/omegam)*(omegac*x1**2-half)*omegap/(lmd**2*sin(theta)**2+omegap**2))\
        * lmd*sin(theta)/(lmd**2*sin(theta)**2+omegap**2)*exp(-lmd*sin(theta)*tau)*cos(omegap*tau)



c_exp2tau=I*omegam/(4*lmd*sin(theta))*(lmd*cos(theta)+Deltam)\
        *(s2+g0*sqrt(2/omegam)*(omegac*x1**2-half)*omegap/(lmd**2*sin(theta)**2+omegap**2))**2*exp(-2*lmd*sin(theta)*tau)


c1=c_tau+c_sin_2omegap+c_cos_2omegap+c_exp_sin+c_exp_cos+c_exp2tau

D1=(lmd**2*sin(theta)**2+omegap**2)
poly_x1=(omegac*x1**2-half)


c_2exp_tau0=I*omegam/(4*lmd*sin(theta))*(lmd*cos(theta)+Deltam)\
    *(s2+g0*sqrt(2/omegam)*poly_x1*omegap/D1)**2


c_cos_2omegap_tau0=I*g0**2/(2*omegap)*lmd*sin(theta)\
    *(lmd**2*sin(theta)**2+omegap**2-omegap*lmd*cos(theta)-omegap*Deltam)/D1**2*poly_x1**2

c_exp_cos_tau0=I*g0*sqrt(2*omegam)*poly_x1*(s2+g0*sqrt(2/omegam)*poly_x1*omegap/D1)*lmd*sin(theta)/D1


#c1 with tau=0
c1_x1_s2_tau0=c_cos_2omegap_tau0+c_exp_cos_tau0+c_2exp_tau0





# z=exp(-c1_x1_s2_tau0)*exp(c1)

# lhs=diff(z,tau)
# rhs=c0*z
#
# tmp=lhs-rhs
# val=tmp.subs([(tau,0.1),(s2,10),(omegap,3),(omegam,2),(omegac,6004),(lmd,2),(g0,10),(theta,12),(Deltam,6),(x1,0.1)])
# pprint(val.evalf())
# pprint(tmp.simplify())

# tmp_func=lambdify((g0,omegam,omegac,omegap,lmd,theta,tau,Deltam,x1,s2),tmp,"numpy")

# valsVec=[20,10,3,2,600,2,10,12,6,0.1]
# print(tmp_func(2,10,3,2,6,2,10,1,6,0.1))
# end: x2 is a function of s2
#################################################

##################################################

X2=symbols("X2",cls=Symbol,real=True)

s2_in_X2=X2*exp(lmd*sin(theta)*tau)\
    -g0*sqrt(2/omegam)*(omegac*x1**2-half)*(lmd*sin(theta)*sin(omegap*tau)-omegap*cos(omegap*tau))/(lmd**2*sin(theta)**2+omegap**2)*exp(lmd*sin(theta)*tau)\
    -g0*sqrt(2/omegam)*(omegac*x1**2-half)*omegap/(lmd**2*sin(theta)**2+omegap**2)


c1_in_X2=c1.subs([(s2,s2_in_X2)])


c1_tau0_in_X2=c1_x1_s2_tau0.subs([(s2,s2_in_X2)])


c_tau_in_X2=c_tau.subs([(s2,s2_in_X2)])
c_cos_2omegap_in_X2=c_cos_2omegap.subs([(s2,s2_in_X2)])
c_sin_2omegap_in_X2=c_sin_2omegap.subs([(s2,s2_in_X2)])
c_exp_cos_in_X2=c_exp_cos.subs([(s2,s2_in_X2)])
c_exp_sin_X2=c_exp_sin.subs([(s2,s2_in_X2)])
c_exp2tau_in_X2=c_exp2tau.subs([(s2,s2_in_X2)])

rho_x1=omegac*x1**2-half
D=lmd**2*sin(theta)**2+omegap**2
mu=lmd*cos(theta)+Deltam

c_tau_short=I*(quarter*omegac+half*Deltam-half*omegac*rho_x1+(2*omegap-mu)/(2*D)*g0**2*rho_x1**2)*tau+half*lmd*sin(theta)*tau


c_cos_2omegap_short=I*half*g0**2*lmd*sin(theta)*(D-mu*omegap)/(omegap*D**2)*rho_x1**2*cos(2*omegap*tau)

c_sin_2omegap_short=I*g0**2*(2*omegap*D+mu*(lmd**2*sin(theta)**2-omegap**2))/(4*omegap*D**2)*rho_x1**2*sin(2*omegap*tau)

c_exp_cos_short=I*g0*sqrt(2*omegam)*lmd*sin(theta)/D*rho_x1*X2*cos(omegap*tau)\
               - I*g0**2*lmd*sin(theta)/D**2*rho_x1**2*(lmd*sin(theta)*sin(2*omegap*tau)-omegap*cos(2*omegap*tau)-omegap)


c_exp_sin_short=I*g0*sqrt(2*omegam)*(mu-omegap)/D*rho_x1*X2*sin(omegap*tau)\
               - I*g0**2*(mu-omegap)/D**2*rho_x1**2*(lmd*sin(theta)-lmd*sin(theta)*cos(2*omegap*tau)-omegap*sin(2*omegap*tau))

c_exp2tau_short=I*omegam*mu/(4*lmd*sin(theta))*(X2-g0*sqrt(2/omegam)*rho_x1*(lmd*sin(theta)*sin(omegap*tau)-omegap*cos(omegap*tau))/D)**2


c_short=c_tau_short+c_cos_2omegap_short+c_sin_2omegap_short+c_exp_cos_short+c_exp_sin_short+c_exp2tau_short

# c_exp2tau_short_expand=I*omegam*mu/(4*lmd*sin(theta))*X2**2\
#               - I*omegam*mu*g0/(2*lmd*sin(theta))*sqrt(2/omegam)*rho_x1*X2*(lmd*sin(theta)*sin(omegap*tau)-omegap*cos(omegap*tau))/D\
#               + I*mu*g0**2/(4*lmd*sin(theta))*rho_x1**2*((omegap**2-lmd**2*sin(theta)**2)*cos(2*omegap*tau)-2*lmd*omegap*sin(theta)*sin(2*omegap*tau))/D**2\
#               +I*mu*g0**2/(4*lmd*D*sin(theta))*rho_x1**2
# tmp=c_exp2tau_short-c_exp2tau_short_expand
# pprint(tmp.simplify())

G_tau=(quarter*omegac+half*Deltam-half*omegac*rho_x1+(2*omegap-mu)/(2*D)*g0**2*rho_x1**2)*tau

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

G=G_tau+G_cos_omegap+G_sin_omegap+G_cos_2omegap+G_sin_2omegap+G_1

func=sin(2*omegap*tau)

c_short_expand=expand(c_short)

P1=quarter*omegac+half*Deltam-half*omegac*rho_x1+(2*omegap-mu)/(2*D)*g0**2*rho_x1**2

F2=g0*sqrt(2*omegam)*(2*lmd**2*sin(theta)**2+omegap*mu)/(2*D*lmd*sin(theta))

F3=g0*sqrt(2*omegam)/D*(half*mu-omegap)

F4=g0**2 * (2*D*lmd**2*sin(theta)**2+lmd**2*mu*omegap*sin(theta)**2+mu*omegap**3)/(4*D**2*lmd*omegap*sin(theta))

F5=g0**2*(2*omegap*D+mu*lmd**2*sin(theta)**2-4*omegap*lmd**2*sin(theta)**2+mu*omegap**2-4*omegap**3)/(4*omegap*D**2)

F6=g0**2*(8*lmd**2*omegap*sin(theta)**2-4*lmd**2*mu*sin(theta)**2+D*mu)/(4*lmd*sin(theta)*D**2)
F7=omegam*mu/(4*lmd*sin(theta))




A=I*P1*tau+I*F2*rho_x1*X2*cos(omegap*tau)+I*F3*rho_x1*X2*sin(omegap*tau)\
    +I*F4*rho_x1**2*cos(2*omegap*tau)+I*F5*rho_x1**2*sin(2*omegap*tau)+I*F6*rho_x1**2+I*F7*X2**2\
    +half*lmd*sin(theta)*tau



B=I*g0**2*lmd*sin(theta)/(2*omegap)*(D-omegap*mu)/D**2*rho_x1**2\
    +I*g0*lmd*sin(theta)/D*sqrt(2*omegam)*rho_x1*X2*exp(lmd*sin(theta)*tau)-I*2*g0**2*lmd**2*sin(theta)**2/D**2*rho_x1**2*sin(omegap*tau)*exp(lmd*sin(theta)*tau)\
    +I*2*g0**2*omegap*lmd*sin(theta)/D**2*rho_x1**2*cos(omegap*tau)*exp(lmd*sin(theta)*tau)\
    +I*omegam/(4*lmd*sin(theta))*mu*X2**2*exp(2*lmd*sin(theta)*tau)+I*mu*g0**2/(4*lmd*sin(theta)*D)*rho_x1**2*exp(2*lmd*sin(theta)*tau)\
    -I*mu*g0/(2*D)*sqrt(2*omegam)*rho_x1*X2*sin(omegap*tau)*exp(2*lmd*sin(theta)*tau)\
    +I*mu*g0*omegap/(2*lmd*sin(theta)*D)*sqrt(2*omegam)*rho_x1*X2*cos(omegap*tau)*exp(2*lmd*sin(theta)*tau)\
    +I*mu*g0**2/(4*lmd*sin(theta)*D**2)*(omegap**2-lmd**2*sin(theta)**2)*rho_x1**2*cos(2*omegap*tau)*exp(2*lmd*sin(theta)*tau)\
    -I*mu*omegap*g0**2/(2*D**2)*rho_x1**2*sin(2*omegap*tau)*exp(2*lmd*sin(theta)*tau)


R1=g0**2*lmd*sin(theta)/(2*omegap)*(D-omegap*mu)/D**2

R2=g0*lmd*sin(theta)/D*sqrt(2*omegam)

R3=-2*g0**2*lmd**2*sin(theta)**2/D**2

R4=2*g0**2*omegap*lmd*sin(theta)/D**2

R5=omegam/(4*lmd*sin(theta))*mu

R6=mu*g0**2/(4*lmd*sin(theta)*D)

R7=-mu*g0/(2*D)*sqrt(2*omegam)

R8=mu*g0*omegap/(2*lmd*sin(theta)*D)*sqrt(2*omegam)

R9=mu*g0**2/(4*lmd*sin(theta)*D**2)*(omegap**2-lmd**2*sin(theta)**2)

R10=-mu*omegap*g0**2/(2*D**2)

B_subs=I*R1*rho_x1**2\
    +I*R2*rho_x1*X2*exp(lmd*sin(theta)*tau)+I*R3*rho_x1**2*sin(omegap*tau)*exp(lmd*sin(theta)*tau)\
    +I*R4*rho_x1**2*cos(omegap*tau)*exp(lmd*sin(theta)*tau)\
    +I*(R5*X2**2+R6*rho_x1**2)*exp(2*lmd*sin(theta)*tau)\
    +I*R7*rho_x1*X2*sin(omegap*tau)*exp(2*lmd*sin(theta)*tau)\
    +I*R8*rho_x1*X2*cos(omegap*tau)*exp(2*lmd*sin(theta)*tau)\
    +I*R9*rho_x1**2*cos(2*omegap*tau)*exp(2*lmd*sin(theta)*tau)\
    +I* R10*rho_x1**2*sin(2*omegap*tau)*exp(2*lmd*sin(theta)*tau)


##################################
# verify pde
z=exp(-B_subs)*exp(A)


c0_in_X2=c0.subs([(x2,X2)])
# pprint(c0_in_X2)

d0=g0*omegac*sqrt(2/omegam)*sin(omegap*tau)*x1**2-half*g0*sqrt(2/omegam)*sin(omegap*tau)-lmd*sin(theta)*X2

lhs=diff(z,tau)+d0*diff(z,X2)

rhs=c0_in_X2*z

tmp=lhs-rhs


val=tmp.subs([(tau,100),(X2,10),(omegap,30),(omegam,2),(omegac,6004),(lmd,2),(g0,10),(theta,12),(Deltam,6),(x1,0.1)])
pprint(val.evalf())
#
##################################
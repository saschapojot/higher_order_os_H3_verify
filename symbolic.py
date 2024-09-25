from sympy import *



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
z=exp(-c1_tau0_in_X2)*exp(c1_in_X2)


c0_in_X2=c0.subs([(x2,X2)])
# pprint(c0_in_X2)

d0=g0*omegac*sqrt(2/omegam)*sin(omegap*tau)*x1**2-half*g0*sqrt(2/omegam)*sin(omegap*tau)-lmd*sin(theta)*X2

lhs=diff(z,tau)+d0*diff(z,X2)

rhs=c0_in_X2*z

tmp=lhs-rhs


val=tmp.subs([(tau,100),(X2,10),(omegap,3),(omegam,2),(omegac,6004),(lmd,2),(g0,10),(theta,12),(Deltam,6),(x1,0.1)])
pprint(val.evalf())
#Based on Coarse-Grain Model for Budding yeast Proposed by
#L. Alberghina1,  E. Martegani, M. Vanoni (University of Milan, Bicocca)
#P. Palumbo (IASI, ROME)


import numpy as np
from scipy import integrate


def yeastode(t,y):

#Input---- the initial data as (n,) numpy array.
#Output---This function gives an (n,1) vector for the RHS of the system of ODE 
# dy/dt = f(y,t).
 
 
 
 
#Allocation of values for caluclating f(y,t) 

  R  = y[0]
  P  = y[1]
  C3F1 = y[2]
  F1 = y[3]
  C3nuc = y[4]
  X1 = y[5]
  X2 = y[6] #Parameter for Timer T2. 
  XB = y[7]
  s=0
  
# Data for growth rate eqn coefficients.


  K0 = 3.0*10**-8
  K1 = 1.0            #1/min
  rho = 1.35*10**-5   #rib/aa
          
  D1 = 4000.0          #min
  D2_0 = 3000
  K2_0 = 542
  
  if s == 0 or s == 1:
  
    D2 = 3000; K2 = 542  #aa/rib/min
    
  elif s==2 :
  
    D2 = 500; K2 = 490
    
  elif s == 3:
  
    D2 = 350; K2 = 430
    
  elif s == 4:
  
    D2= 240;K2= 350
    
  elif s == 5:
  
    D2 ==220 ; K2= 345
    
  else:
  
    D2 = 220 ; K2 = 340
    
# Initialize for paramater X2 (Refer coeff of P(t) ode, 
# we need X2 for this calculation of coeff)

  beta2  = 1.0
  gamma2 = 0.01
  T2     = 15.0
  n2     = 100
  X2_thr = (beta2/gamma2)*(1.0 - np.exp(-gamma2*T2)) # Threshold of for X2.

  
  K2_f = (K2_0 - K2)*((X2/X2_thr)**n2/(1 + (X2/X2_thr)**n2)) +  K2
  D2_f = (D2_0 - D2)*((X2/X2_thr)**n2/(1 + (X2/X2_thr)**n2)) +  D2

#Data for parameter X1

  gamma1 = 0.01 #[min^-1]
  beta1 = 1.0 #[au/min]
  n1 = 100
  lambda_1 = rho*K2 - 1/D2 
  
  if s == 0:
    T1 = (1/lambda_1)*np.log(2.50/1.68) - T2 #for new born daughtercell.
  else:
    T1 = 1  #for parent cells s=1,2,...
  X1_thr = (beta1/gamma1)*(1.0 - np.exp(-gamma1*T1))
  
  
#Data for XB paramater

  gamma_B = 0.01    #[min^-1]
  beta_B = 1.0      #[au/min]
  n_B = 100
  T_B = 74          #[min]
  T_B_star = 0.8*T_B
  XB_thr = (beta_B/gamma_B)*(1.0 - np.exp(-gamma_B*T_B))
  XB_star = (beta_B/gamma_B)*(1.0 - np.exp(-gamma_B*T_B_star))

#Data for SIZER 

  K_on  = 1.6284*10**-15 #[(mol/L)^-1/min]
  K_off = 25 #[min^-1]
  h = .07
  H = 7.09*10**23 # [aa/L]
  V_tot = P/H
  V_nuc = h*V_tot
  n_f = 10
  K_nc = 0.6 #[min^-1]
  theta = 3.02*10**-8 #[mol/aa]
  C3_tot = theta*P
  
  if XB <= XB_star :
    eta = 1.0 #[min^-1]
  else :
    eta = 0.0
    
  if XB <= XB_star :
    K_cn = 1.5 #[min^-1]
  else :
    K_cn = 5*10**-4
  
  if XB <= XB_star :
    alpha = 1.0 
  else :
    alpha = 2.0
    
  C3_cyt = C3_tot - alpha*C3nuc + C3F1

  
#output function


  n=len(y)
  dydt=np.zeros((n,1))
  if rho*P > R:
    dydt[0] = K0*P + K1*(rho*P - R) - R/D1
  else:
    dydt[0] = K0*P - R/D1

  dydt[1] = K2_f*R - P/D2_f
  dydt[2] = (K_on/V_nuc)*(C3nuc*F1) - K_off*C3F1
  dydt[3]= -(K_on/V_nuc)*(C3nuc*F1) + K_off*C3F1 - eta*((C3nuc/C3F1)**n_f/(1+(C3nuc/C3F1)**n_f))*F1
  dydt[4]= -(K_on/V_nuc)*(C3nuc*F1) + K_off*C3F1 + K_cn*C3_cyt - K_nc*C3nuc
  dydt[5]= -gamma1*X1 + beta1*(1/(1+(C3F1/C3nuc)**n1))
  dydt[6] = -gamma2*X2 + beta2*(1/(1+(X1_thr/X1)**n2))
  dydt[7]= -gamma_B*XB + beta_B*(1/(1+(X2_thr/X2)**n_B))
  return dydt
  
  
def main():
  #Choice of integrator.
  
  
  r = integrate.ode(yeastode).set_integrator('vode',method='bdf')
  
  
  #Time steps
  
  t_s = 0.0
  t_f = 135.0
  delta_t = 0.1
  
  
  num_steps = np.floor((t_f - t_s)/delta_t)+1
  
  #Initialize ODE
  
  R_t_zero = 2.25*10**5
  P_t_zero = 1.68*10**10
  C3F1_t_zero = 0.0001
  F1_t_zero = 238.0
  C3nuc_t_zero = 9.0
  X1_t_zero = 0.0001
  X2_t_zero = 0.0001
  XB_t_zero = 0.0001
  
  r.set_initial_value([R_t_zero,P_t_zero,C3F1_t_zero,F1_t_zero,C3nuc_t_zero,X1_t_zero,X2_t_zero,XB_t_zero],t_s)
  
  #Array to store data
  
  t = np.zeros((num_steps,1))
  rib=np.zeros((num_steps,1))
  pro=np.zeros((num_steps,1))
  c3f1 = np.zeros((num_steps,1))
  f1 = np.zeros((num_steps,1))
  c3nuc = np.zeros((num_steps,1))
  x1 = np.zeros((num_steps,1))
  x2 = np.zeros((num_steps,1))
  xb = np.zeros((num_steps,1))
  
  t[0]=t_s
  rib[0]=R_t_zero
  pro[0]=P_t_zero
  c3f1[0] = C3F1_t_zero
  f1[0] = F1_t_zero
  c3nuc[0] = C3nuc_t_zero
  x1[0] = X1_t_zero
  x2[0] = X2_t_zero
  xb[0] = XB_t_zero
  
  #Running the integration.
  
  z = 1
  while r.successful() and z < num_steps:
  
    r.integrate(r.t + delta_t)
    
    t[z]=r.t
    rib[z] = r.y[0]
    pro[z]=r.y[1]
    c3f1[z] = r.y[2]
    f1[z] = r.y[3]
    c3nuc[z] = r.y[4]
    x1[z] = r.y[5]
    x2[z] = r.y[6]
    xb[z] = r.y[7]
    z += 1

#Plotting of the parameters.
    
  import pylab as py
  py.subplot(211)
  py.plot(t,rib,'r-')
  py.xlabel('Time')
  py.ylabel('Ribosomes #rib ')
  py.grid('on')
  py.legend()
  py.title('Growth , s=0')

  py.subplot(212)
  py.plot(t,pro,'g-')
  py.xlabel('time')
  py.ylabel('Protein #aa')
  py.legend()
  py.grid('on')
  
  
  '''py.subplot(321)
  py.plot(t,c3f1,label='Cln3Far1')
  py.xlabel('time')
  py.ylabel('Cln3Far1')
  py.legend()
  py.grid('on')
  
  py.subplot(322)
  py.plot(t,f1,label='Far1')
  py.xlabel('Time')
  py.ylabel('Far1')
  py.legend()
  py.grid('on')
  
  py.subplot(323)
  py.plot(t,c3nuc,label='Cln3nuc')
  py.legend()
  py.xlabel('Time')
  py.grid('on')
  
  py.subplot(324)
  py.plot(t,x1,label='X1')
  py.xlabel('Time')
  py.ylabel(' X1')
  py.xlim([0,40])
  py.legend()
  py.grid('on')

  py.subplot(325)
  py.plot(t,x2,label='X2')
  py.xlabel('Time')
  py.ylabel(' X2')
  py.legend()
  py.xlim([40,55])
  py.grid('on')
  
  
  py.subplot(326)
  py.plot(t,xb,label='XB')
  py.xlabel('Time')
  py.ylabel(' X_B')
  py.legend()
  py.grid('on')'''
  
  py.show()
    
  


#Boiler plate syntax
if __name__ == '__main__':
  main()



# coding: utf-8

# A simplified model similar to [O'Leary et al. (2014)](http://alexhwilliams.info/pubs/OLeary_etal_2014.pdf)
# -------------
# 
# This is a model of activity-dependent conductance regulation, based on the time-averaged intracellular calcium concentration. While this can (and has) be used to regulate voltage-dependent conductances with Hodgkin-Huxley type channels, a simpler conductance-based model with only leak currents suffices to illustrate the basic behavior of the model. This is because the model integrates the average calcium concentration over a fairly large time window, and is therefore insensitive to calcium transients produced in realistic models (but see [Liu et al., 1998](http://neurotheory.columbia.edu/~larry/LiuJNeurosci98.pdf)).
# 
# The code below numerically simulates the following set of equations.
# $$ \frac{dV}{dt} = \sum_i g_i (E_i - V)~~~~~~~~~~~~~~(1)$$
# $$ c(V) = \frac{1}{1+\exp(-V/20)}~~~~~~~~~~(2)$$
# $$ ~~~~ \frac{dg_i}{dt} = \frac{u_i - g_i}{\tau_g}~~~~~~~~~~~~~~~~~~~~~~~(3)$$
# $$ ~~~~~~ \frac{du_i}{dt} = \alpha_i(T) - \beta_i u_i~~~~~~~~~~~~~~~~~~(4)$$
# $$ ~~~~~~~ \frac{dT}{dt} = \alpha_T(c) - \beta_T~~~~~~~~~~~~~~~~~~~~~~(5)$$
# 
# Equation (1) is the membrane equation. Equation (2) specifies the intracellular calcium level, which is approximated as a monotonic function of membrane potential. Equation (3) specifies the change in each leak conductance, $g_i$. Note that $g_i \rightarrow u_i$ over time, where $u_i$ is the 'mRNA concentration' for that ion channel type. We assume this process occurs for all ion channels with the same time constant, $\tau_g$. It is possible to instead solve this equation with channel-specific insertion and removal rate constants:
# 
# $$ ~~~~ \frac{dg_i}{dt} = \alpha_{g,i} u_i - \beta_{g,i} g_i~~~~~~~~~~~~~~~~~~~~~~~(3b)$$
# 
# However, the simplified form of equation (3) cuts down the parameters in the model without modifying its overall behavior. Equation (4) specifies the change in mRNA concentration over time for each channel type, where $\alpha_i(T)$ is the production rate and $\beta_i$ proportional to the degradation rate. The production of each mRNA is determined by the concentration of a 'master' transcription factor $T$, whose dynamics are provided by equation (5); $\alpha_T$ and $\beta_T$ are constants. Assuming constant degradation is necessary for conductance regulation to implement integral control (O'Leary et al., 2014). It is reasonable to model $\alpha_i(T)$ in equation (4) as well as $\alpha_T(c)$ in equation (5) using scaled Hill Equations:
# 
# $$ ~~~~ \alpha_i(T) = \alpha_{min} + (\alpha_{max} - \alpha_{min}) \left ( \frac{T^n}{T^n + K_d} \right ) ~~~~~~~~~~~~~~~~~~~~~~~(6a)$$
# $$ ~~~ \alpha_i(T) = \alpha_{min} + (\alpha_{max} - \alpha_{min}) \left ( 1-\frac{T^n}{T^n + K_d} \right ) ~~~~~~~~~~~~~~~~~~~~(6b)$$
# 
# Equation (6a) describes a scenario where the binding of $T$ induces transcription, while (6b) describes a scenario where $T$ suppresses transcription. In both cases $\alpha_{max}$ and $\alpha_{min}$ are the soft bounds on the transcription rate. Note that these equations place bounds on $u_i$, which precludes the system from performing 'perfect' integral control because the integrated error signal can saturate.
# 
# Note that the 'target', 'set point', or steady-state calcium concentration, $c_{ss}$, is found by solving:
# $$ \alpha_T (c_{ss}) = \beta_T $$
# When $c = c_{ss}$, the system is in equilibrium, since $dT/dt = 0$.

# In[1]:

## This cell determines all the parameters in the model
from __future__ import division
import numpy as np
import pylab as plt
get_ipython().magic(u'matplotlib inline')

# NOTE: lambda is the syntax that produces an anonymous functoin
#  --> http://www.secnetix.de/olli/Python/lambda_functions.hawk

# Parameters for equations 1-3
E = [-80, -10, 50]  # list of reversal potentials each leak current
tau_g = 10.0        # arbitrary time units
c = lambda V: 1/(1+(np.exp(-V/20))) # function for calcium, e.g. c(0) yields 0.5

# Equation 4, with alpha given by equations 6a & 6b
beta = [1,1,1]
alpha = []
alpha.append(lambda T: 0.0 + 2*(1-(T**2/(T**2+1)))) 
alpha.append(lambda T: 0.0 + 4*(T/(T+4)))
alpha.append(lambda T: 0.0 + 2*(T**2/(T**2+1)))

# Equation 5, steady state at c = beta_T / alpha_T
alpha_T = lambda ca: 0.1*(1-(ca**2/(ca**2+0.2)))
beta_T = 0.05


# In[2]:

## This cell plots some useful information about the parameters

# Get a nice colormap for plotting
# --> brewer2mpl library is here: https://github.com/jiffyclub/brewer2mpl
import brewer2mpl
bmap = brewer2mpl.get_map('Set1', 'qualitative', len(E))
colors = bmap.mpl_colors

v_axis = np.linspace(min(E),max(E),100)
T_axis = np.linspace(0,5,100)

# Calcium as a function of membrane potential
plt.figure()
plt.plot(v_axis,c(v_axis),lw=2)
plt.xlabel('membrane potential (mV)',fontweight='bold')
plt.ylabel('Ca conc. (a.u.)',fontweight='bold')
plt.show()

# Transcription rate as a function of master transcription factor for each channel type
# --> legend labels each channel type by the reversal potential
plt.figure()
for (i,a) in enumerate(alpha):
    plt.plot(T_axis,a(T_axis),'-',c=colors[i],lw=2,label=str(E[i]))
plt.xlabel('conc. of T',fontweight='bold')
plt.ylabel('transcription rate',fontweight='bold')
plt.legend(loc='right')
plt.show()

# Transcription factor activation rate as a function of calcium
plt.figure()
plt.plot(c(v_axis),alpha_T(c(v_axis))-beta_T,'-k',lw=2)
plt.plot(c(v_axis),np.zeros(len(v_axis)),'--r',lw=2)
plt.xlabel('Ca conc. (a.u.)',fontweight='bold')
plt.ylabel('dT/dt',fontweight='bold')
plt.show()


# **Notes on the above plots**
# 
# * The target calcium concentration or homeostatic set point for the system occurs when $dT/dt$ is equal to zero. This occurs at the intersection of the black line and red dashed line in the above plot (about 0.5 in this case). Note that this fixed point is stable because the $d^2T/dt^2 < 0$ at that point (i.e. the slope of the black line is negative).
# 
# * If $dT/dt$ is nearly linear, and doesn't saturate, then the system approximately implements integral control.
# 
# **Two final notes on this implementation**
# 
# 1) Rather than numerically integrate equation (1) explicitly, we set the membrane potential equal to the weighted sum of the reversal potentials. This is Ohm's Law:
# $$V = V_{ss} = \frac{\sum{g_i E_i}}{\sum{g_i}} $$
# 
# 2) The model allows $T$ to go negative. However, I fix the input to $\alpha_i(\cdot)$ to be $\max(0,T)$.

# In[3]:

## This cell does the simulation
from scipy.integrate import odeint # import ODE solver

n = len(E) # the number of ion channel types

# steady-state membrane potential
def v(g):
    return np.sum(g*E)/np.sum(g)

# a rectifier function to prevent T<0
def fix(x):
    if x < 0: return 0
    else: return x

# function that returns dy/dt given vector y as input
def model(y,t0):
    # y[0:n] are the conductances, g
    # y[n:2*n] are the mRNA concentrations, u
    # y[-1], the last element, is T
    
    dydt = []  # list of derivatives, vars in same order as y
    
    Vnow = v(y[0:n]) # current membrane potential
    c_now = c(Vnow)  # current calcium concentration
    T = fix(y[-1])   # current value of T (not negative)
    
    # Change in conductance
    for i in range(n):
        dydt.append((y[i+n]-y[i])/tau_g) # Equation (3)
    
    # Change in mRNA concentration
    for (i,a) in enumerate(alpha):
        dydt.append(a(T)-beta[i]*y[i+n]) # Equation (4)
    
    # Change in T
    dydt.append(alpha_T(c_now) - beta_T)
    
    # Return calculated derivatives
    return dydt


# initial conductances
g0 = np.array([2,0.1,0.2]) + np.random.uniform(-0.1,0.1,3)
u0 = g0          # initial mRNA
T0 = [1.0]       # initial transcription factor

# initial state
y0 = np.concatenate((g0,u0,T0))
t = np.linspace(0,100)

# Numerically solve the equations given by the function 'model'
y = odeint(model, y0, t)

g = y[:,0:n]    # conductances over time
u = y[:,n:2*n]  # mRNA over time
T = y[:,-1]     # T over time

# membrane potential over time
v = np.dot(g,E)/np.sum(g,axis=1)


# In[4]:

## More plots

# Conductances over time
plt.figure()
for i in range(3):
    plt.plot(t,g[:,i],c=colors[i],lw=2,label=str(E[i]))
plt.legend()
plt.ylabel('conductances',fontweight='bold')
plt.show()

# membrane potential over time
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(t, v, 'k-',lw=2)
ax2.plot(t, c(v), 'b-',lw=2)
ax1.set_xlabel('time (a.u.)',fontweight='bold')
ax1.set_ylabel('membrane potential (mV)', color='k',fontweight='bold')
ax2.set_ylabel('Ca conc. (a.u.)', color='b',fontweight='bold')
plt.show()


# **Notes on the above plots:**
# 
# * The initial condition has a large conductance for the "red channel type" (middle plot, above) with a reversal potential of -80mV. This results in a calcium concentration below the set point. Over time, the conductance of the red channel type decreases as the other two increase. Eventually the membrane potential and the calcium concentration increase and stabilize at their target values (bottom plot).

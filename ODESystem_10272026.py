#!/usr/bin/env python
# coding: utf-8

# #### Snyder ODE System

# In[3]:


def SnyderODE(t, N,                                         #timesteps
              rC, rM, rT, αCT, αMT, αMC, μC,                #benthos
              sB, hB, zB, μB,                               #browsers
              sG, hG, zG, μG,                               #grazers
              sS, hS, zS, μS):                              #scrapers   

    
    C, M, T, B, G, S = N

    
    #Hollings_B = (hB*B*M)/((hB*ηB*M)+1)
    #Hollings_G = (hG*G*T)/((hG*ηG*T)+1)
    #Hollings_S = (hS*S*T)/((hS*ηS*T)+1)

    Blackwood_B = ((1-zB)*C) + (zB*M)
    Blackwood_G = ((1-zG)*C) + (zG*T)
    Blackwood_S = ((1-zS)*C) + (zS*T)
    
    K_B = (1-(B/(0.001 + ((0.0564 - 0.001)*Blackwood_B)))) 
    K_G = (1-(G/(0.001 + ((0.0564 - 0.001)*Blackwood_G))))
    K_S = (1-(S/(0.001 + ((0.0564 - 0.001)*Blackwood_S))))
    
    
    dC = rC*C*(1-C-M-T) + αCT*rC*C*T - αMC*rM*M*C - μC*C
    dM = rM*M*(1-C-M-T) + αMT*rM*M*T + αMC*rM*M*C - hB*B*M #Hollings_B    
    dT = rT*T*(1-C-M-T) - αCT*rC*C*T - αMT*rM*M*T - hG*G*T - hS*S*T #Hollings_G - Hollings_S
    
    dB = (sB*B*K_B) - μB*B
    dG = (sG*G*K_G) - μG*G
    dS = (sS*S*K_S) - μS*S 

      
    return dC, dM, dT, dB, dG, dS


# In[ ]:





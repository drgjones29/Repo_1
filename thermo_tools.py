from math import pow, log,exp
import numpy as np
from parse_chem import * 

#define some constants used here
R = 8.31450  # JK-1mol
h = 6.63e-34   # Js
kb= 1.38e-23   # JK-1
c = 3.00e+10   # cms-1
L = 6.022e+23  # 
hc_over_k = h*c/kb 
kbeV = 8.617269e-5 # eVK-1

check=open('check.out','w')       
    
def generate_species_list(reactions):
    
    species_list={}
    l=len(reactions)
    
    for part in['lhs','rhs','TS']:
        for i in range(0,l,1):    
	    
	        try:
                    spec_dum = reactions.get('r'+str(i))[part]
		         
	            #print 'spec dum', spec_dum
	            for j in spec_dum :
	
	                #print j
	    
	                if j not in species_list.itervalues():
	        
		            species_list[len(species_list)]=j 
                except KeyError:
		    print >> check, "No Key found for part "+ part +" in reaction "+ str(i) 

    return species_list
    
def get_G_species(species,datasurf=None,datagas=None,Temp=298.,S_gas='shom',S_surf=None, \
                      ZPE_gas=None,ZPE_surf=None,H_gas=None,U_surf=None,P_work=None,over_pot=None, approx=None):
    
    if approx == 'edft':
        print 'Returning plain dft (with an ad hoc adjsutments you make'
        S_gas=None ; S_surf=None ; ZPE_gas=None ; ZPE_surf=None ; H_gas=None ; U_surf=None ; P_work=None ; over_pot=None
    elif approx == 'g_ent_only':
        print 'Using the shomate fitted entropy only for the gas_phase'
        S_gas='shom' ; S_surf=None ; ZPE_gas=None ; ZPE_surf=None ; H_gas=None ; U_surf=None ; P_work=None ; over_pot=None
    elif approx == 'harm_st':
        print 'Using harmonic theory of a classical ideal gas for surf'
        S_gas='shom' ; S_surf='vib' ; ZPE_gas='zpe' ; ZPE_surf='vib' ; H_gas='shom' ; U_surf='vib' ; P_work=None ; over_pot=None
    elif approx == 'g_trans':
        print 'Using the translation entropy only for the gas_phase'
        S_gas='trans' ; S_surf=None ; ZPE_gas=None ; ZPE_surf=None ; H_gas=None ; U_surf=None ; P_work=None ; over_pot=None 
    else:
        print 'Using default settings to calcualte G, is this what you are expecting?'    
    #at the moment specific choices for each thermodynamic contribution need to be made (otherwise they are ignored and E is returned
    #At a later date coudl introduce default approximations too e.g. approx1=(blardy blah), approx2=another set of values
    
    species_G=np.zeros(len(species))
    
    for i in range(0,len(species),1):
                   
       phase=species.get(i)[-1] 
       print 'Running for SPECIES ** :' , species.get(i)
       print 'T = ', Temp
       if phase == '*':
           #print 'Species get', species.get(i)
      	   spec = datasurf[0].get(datasurf[1]).get(species.get(i))
	   
	   #Epot
	   	   
	   if datasurf[0] == 'scaling' :
	      print >> check, " trying to implement differently "
	      
	   else:
	      	       
	       E_pot = spec['Edft']
	       print "getting Edft!", E_pot
	       
	   try:
	       Fudge = spec['Fudge']	
	       if Fudge is not None:
	           print "adding afudge factor ! warning!"
	           E_pot=E_pot+ Fudge   	
		   print "E_pot", E_pot
	       else:
	           print 'no ad hoc corrction'	          
	   except KeyError:
	       print 'no fudge factor found' 
	       
	   #if Fudge is not None:      
	   #    E_pot=E_pot+ Fudge	   
		         
	   #ZPE
	   if ZPE_surf == None : 
	       zpe=0.0
	       print 'no zpe added'
	   elif ZPE_surf == 'vib' :
	       try: 	       	       
	           zpe=calc_zpe_vib(spec['vib'])
		   print 'vibrational zpe:',zpe
	       except KeyError:
	           print >> check, "You don't appear to have any vibrations of species "+species.get(i)
		   print >> check, "zpe set to 0.0"
	           print >> check, "I would try the zero entropy and zpe surface model" 
		   print >> check, "or do mmore DFT calcualtions" 
	   elif ZPE_surf == 'zpe' :      
	       print >> check, "this isn't implemented yet, zpe set to 0.0"
	       print >> check, "I would try the zero entropy and zpe surface model"
	       print >> check, "or do more DFT calculations"
	       zpe=0.0
	       
	   #deltaU
	   if U_surf == None :
	       U= 0.0
	       print 'U surf :', U
	   elif U_surf == 'vib' :   
	       U=calc_u_surf_vib(Temp,spec['vib'])
	       print 'U vib :', U
	   #entropy  
	   if S_surf == None :
	       S = 0.0
	       print 'S surf :', S
	   elif S_surf == 'vib':
	       S=calc_s_surf_vib(Temp,spec['vib'])
               print 'S vib :', S
	   #over potential          
	   if over_pot == None:
	       U_app = 0.0
	       print 'U_app :', U_app
	   else:
	       print >> check, 'applied potential is not implemented yet'
	   #fudge factors determiend from database or dictionary 
	   #you should be able to add a fudge - but need to work out the syntax
           #for now fudge factors will just be modified in the database but a better way may be forthcoming......
	   print >> check, 'calcualting G :'
	   print >> check, species.get(i),E_pot,zpe,U,(Temp*S),U_app 
	   print 'calcualting G :'
	   print species.get(i),E_pot,zpe,U,(Temp*S),U_app  
	   print E_pot+zpe+U-(Temp*S)+U_app
	   species_G[i]=E_pot+zpe+U-(Temp*S)+U_app 
	           
       elif phase == 'g':
           #print 'species.get(i)', species.get(i)
           spec = datagas[0].get(datagas[1]).get(species.get(i))
           print >> check, 'spec',spec          	  
	   #Epot	       	       
	   E_pot = spec['Edft']
	   #print 'test', species.get(i),E_pot
	   
	   try:
	       Fudge = spec['Fudge']
	       if Fudge is not None:	           
	       	   print "adding afudge factor ! warning!"
	           E_pot=E_pot+ Fudge   	
		   print "E_pot", E_pot		   	       
	   except KeyError:
	       print >> check,  'no fudge factor found' 
	       
	   #if Fudge is not None:      
	   #    E_pot=E_pot+ Fudge
	                
           #print 'got here 1 - for debugging'
	   #ZPE
	   if ZPE_gas == None : 
	       zpe=0.0
	       print 'no zpe added'
	   elif ZPE_gas == 'vib' :
	       try: 	       
	           zpe=calc_zpe_vib(spec['vib'])
		   print 'vibrational zpe:',zpe
	       except KeyError: 
	           print >> check, "You don;t appeat to have any vibrations for species "+species.get(i)
		   print >> check, "zpe set to 0.0"
	           print >> check, "I would try the zero entropy and zpe surface model"
		   zpe=0	   
	   elif ZPE_gas == 'zpe' :
	       
	       try:      	       
	           zpe=spec['ZPE']
		   print 'pre_calculated zpe:',zpe
	       except KeyError:
	           print >> check, "you don't appear to have a ZPE for species "+species.get(i)
		   print >> check, "zpe set to 0.0"
	           print >> check, "I would try a different approximation or do some more DFT calcualtions"
		   zpe=0
	   #deltaH
	   if H_gas==None:
	       H=0.0
	       print 'H gas corr.:',H
	   elif H_gas == 'vib': 
	       try:           
	           H=calc_h_gas_vib(Temp,spec['vib'])
		   print 'H vib:', H
	       except KeyError:  
	           print >> check, "You don;t appeat to have any vibrations for species "+species.get(i)
		   print >> check, "deltaH set to 0.0, I would try shomate"		   
		   H=0.0 
	   elif H_gas == 'shom':
	       try:
	           H=calc_h_gas_shom(Temp,spec['shomate'])
		   print ' H shomate:', H
	       except KeyError:
	           print >> check, "you appear to missing the shomate parameters"
		   print >> check, "delta H is being set to zero for species " + species.get(i)
		   print >> check, "try looking them up , or fitting a polynomial to tabulated data"
		   H=0.0  
	   if H_gas is not None:
	       try:
	           H_273=spec['dH_273'] 
		   print 'H 273 :', H_273/96.485
		   H=H+(H_273/96.485)
	       except KeyError:
	          print >> check, "you appear to missing the H_273 value"
		  print >> check, "delta H is being set to zero for species " + species.get(i)
		  print >> check, "try looking them up , or fitting a polynomial to tabulated data"
		   	       	   
	   #entropy  
	   if S_gas == None :
	       S = 0.0
	       print 'S gas:', S
	   elif S_gas == 'vib':
	       try:
	           S=calc_s_gas_vib(Temp,spec['vib'])
		   print 'S vib:', S
	       except KeyError:  
	           print >> check, "You don;t appear to have any vibrations for species "+species.get(i)
		   print >> check, "deltaS set to 0.0, I would try shomate"
		   S=0.0
           elif S_gas == 'shom':
	       try:
	           S=calc_s_gas_shom(Temp,spec['shomate'])
		   print 'S shomate:', S
	       except KeyError:
	           print >> check, "you appear to missing the shomate parameters"
		   print >> check, "delta S is being set to zero for species " + species.get(i)
		   print >> check, "try looking them up , or fitting a polynomial to tabulated data"
		   S=0.0
           elif S_gas == 'trans':
	       #print 'calcualting S'
	       try:
	           
	           S=calc_s_gas_trans(Temp,spec['M'])
		   print 'S trans:', S		  
	       except KeyError:
	           print >> check, "you appear to missing the molecular mass"
		   print >> check, "delta S is being set to zero for species " + species.get(i)
		   print >> check, "There's a periodic table behind you ;)" 	
           	       	   
	   #P_work
	   if P_work == None:
	       p_w = 0.0
	       print 'pressure work is:', p_w
	   else: 
	       try:
	           p_w=calc_p_work(Temp,P_work[species.get(i)])
		   print 'kTln(p/pstan):', p_w
	       except KeyError:
	           print >> check, "you appear to missing something todo with the pressure"
		   print >> check, "Check what you have entered for " + species.get(i)		   
		   p_w=0.0
	   #overpotential (or applied potential)    
	   if over_pot == None:
	       U_app = 0.0
	       print 'Applied pot :' , U_app
	   else:
	       print >> check, 'applied potential is not implemented yet' 
	   print >> check, species.get(i),E_pot, zpe, H, (Temp*S), p_w, U_app, E_pot+zpe+H-(Temp*S)+p_w+U_app 
	   
	   
	   print 'calculating G:'
	   print  species.get(i),E_pot, zpe, H, (Temp*S), p_w, U_app  
	   print  E_pot+zpe+H-(Temp*S)+p_w+U_app    
           species_G[i]=E_pot+zpe+H-(Temp*S)+p_w+U_app
	   
           
    return species_G

def get_H_species(species,datasurf=None,datagas=None,Temp=298., H_gas='shom',U_surf=None):
    #at the moment specific choices for each thermodynamic contribution need to be made (otherwise they are ignored and E is returned
    #At a later date coudl introduce default approximations too e.g. approx1=(blardy blah), approx2=another set of values
  species_H=np.zeros(len(species))
  
  for i in range(0,len(species),1):
    print i  
    print'**************************'
    print 'calcualting H for species', species.get(i) 
    phase=species.get(i)[-1]
    if phase == '*':
           print 'Species get', species.get(i)
      	   spec = datasurf[0].get(datasurf[1]).get(species.get(i))
	   
	   #Epot
	   	   
	   if datasurf[0] == 'scaling' :
	      print >> check, " trying to implement differently "
	      
	   else:
	       	       
	       E_pot = spec['Edft']
	       print "getting Edft!", E_pot
	       
	   try:
	       Fudge = spec['Fudge']	
	       if Fudge is not None:
	           print "adding afudge factor ! warning!"
	           E_pot=E_pot+ Fudge   	
		   print "E_pot", E_pot
	       else:
	           print 'no ad hoc corrction'	          
	   except KeyError:
	       print 'no fudge factor found'  
	   	         	  
	       
	   #deltaU
	   if U_surf == None :
	       U= 0.0
	       print 'U surf :', U
	   elif U_surf == 'vib' :   
	       U=calc_u_surf_vib(Temp,spec['vib'])
	       print 'U vib :', U    
	   species_H[i]=E_pot+U 
	           
    elif phase == 'g':
           
           spec = datagas[0].get(datagas[1]).get(species.get(i))
           print 'spec',spec          	  
	   #Epot	       	       
	   E_pot = spec['Edft']
	   try:
	       Fudge = spec['Fudge']
	       if Fudge is not None:	           
	       	   print "adding afudge factor ! warning!"
	           E_pot=E_pot+ Fudge   	
		   print "E_pot", E_pot		   	       
	   except KeyError:
	       print >> check,  'no fudge factor found'   
	            
        
	   #deltaH
	   if H_gas==None:
	       H=0.0
	       print 'H gas corr.:',H
	   elif H_gas == 'vib': 
	       try:           
	           H=calc_h_gas_vib(Temp,spec['vib'])
		   print 'H vib:', H
	       except KeyError:  
	           print >> check, "You don;t appeat to have any vibrations for species "+species.get(i)
		   print >> check, "deltaH set to 0.0, I would try shomate"		   
		   H=0.0 
	   elif H_gas == 'shom':
	       try:
	           H=calc_h_gas_shom(Temp,spec['shomate'])
		   print ' H shomate:', H
	       except KeyError:
	           print >> check, "you appear to missing the shomate parameters"
		   print >> check, "delta H is being set to zero for species " + species.get(i)
		   print >> check, "try looking them up , or fitting a polynomial to tabulated data"
		   H=0.0  
	   if H_gas is not None:
	       try:
	           H_273=spec['dH_273'] 
		   print 'H 273 :', H_273/96.485
		   H=H+(H_273/96.485)
	       except KeyError:
	          print >> check, "you appear to missing the H_273 value"
		  print >> check, "delta H is being set to zero for species " + species.get(i)
		  print >> check, "try looking them up , or fitting a polynomial to tabulated data"
		   	       	   	   	       
           species_H[i]=E_pot+H
           
  return species_H

def get_S_species(species,datasurf=None,datagas=None,Temp=298., S_gas='shom',S_surf=None):
    #at the moment specific choices for each thermodynamic contribution need to be made (otherwise they are ignored and E is returned
    #At a later date coudl introduce default approximations too e.g. approx1=(blardy blah), approx2=another set of values
    
    species_S=np.zeros(len(species))
    
    for i in range(0,len(species),1):
                   
       phase=species.get(i)[-1] 
       print >> check, 'Running for SPECIES ** :' , species.get(i)
       
       if phase == '*':
      	   spec = datasurf[0].get(datasurf[1]).get(species.get(i))
	   
	   
	       
	  
	   #entropy  
	   if S_surf == None :
	       S = 0.0
	       print 'S surf :', S
	   elif S_surf == 'vib':
	       S=calc_s_surf_vib(Temp,spec['vib'])
               print 'S vib :', S
	   	    	           
       elif phase == 'g':
           
           spec = datagas[0].get(datagas[1]).get(species.get(i))
           print 'spec',spec          	  
	   	  	   
	   #entropy  
	   if S_gas == None :
	       S = 0.0
	       print 'S gas:', S
	   elif S_gas == 'vib':
	       try:
	           S=calc_s_gas_vib(Temp,spec['vib'])
		   print 'S vib:', S
	       except KeyError:  
	           print >> check, "You don;t appear to have any vibrations for species "+species.get(i)
		   print >> check, "deltaS set to 0.0, I would try shomate"
		   S=0.0
           elif S_gas == 'shom':
	       try:
	           S=calc_s_gas_shom(Temp,spec['shomate'])
		   print 'S shomate:', S
	       except KeyError:
	           print >> check, "you appear to missing the shomate parameters"
		   print >> check, "delta S is being set to zero for species " + species.get(i)
		   print >> check, "try looking them up , or fitting a polynomial to tabulated data"
		   S=0.0
           elif S_gas == 'trans':
	       #print 'calcualting S'
	       try:
	           
	           S=calc_s_gas_trans(Temp,spec['M'])
		   print 'S trans:', S		  
	       except KeyError:
	           print >> check, "you appear to missing the molecular mass"
		   print >> check, "delta S is being set to zero for species " + species.get(i)
		   print >> check, "There's a periodic table behind you ;)"	   
	   	         
           species_S[i]=(S)
           
    return species_S 

def get_ZPE_species(species,datasurf=None,datagas=None,Temp=298., ZPE_gas='ZPE',ZPE_surf=None):
    #at the moment specific choices for each thermodynamic contribution need to be made (otherwise they are ignored and E is returned
    #At a later date coudl introduce default approximations too e.g. approx1=(blardy blah), approx2=another set of values
    
    species_ZPE=np.zeros(len(species))
    
    for i in range(0,len(species),1):
                   
       phase=species.get(i)[-1] 
       print >> check, 'Running for SPECIES ** :' , species.get(i)
       
       if phase == '*':
      	   spec = datasurf[0].get(datasurf[1]).get(species.get(i))
	   	   		         
	   #ZPE
	   if ZPE_surf == None : 
	       zpe=0.0
	       print 'no zpe added'
	   elif ZPE_surf == 'vib' :
	       try: 	       	       
	           zpe=calc_zpe_vib(spec['vib'])
		   print 'vibrational zpe:',zpe
	       except KeyError:
	           print >> check, "You don't appear to have any vibrations of species "+species.get(i)
		   print >> check, "zpe set to 0.0"
	           print >> check, "I would try the zero entropy and zpe surface model" 
		   print >> check, "or do mmore DFT calcualtions" 
	   elif ZPE_surf == 'zpe' :      
	       print >> check, "this isn't implemented yet, zpe set to 0.0"
	       print >> check, "I would try the zero entropy and zpe surface model"
	       print >> check, "or do more DFT calculations"
	       zpe=0.0
	       	  	      
	   species_ZPE[i]=zpe 
	           
       elif phase == 'g':
           
           spec = datagas[0].get(datagas[1]).get(species.get(i))
           print 'spec',spec          	  
	   
	   #ZPE
	   if ZPE_gas == None : 
	       zpe=0.0
	       print 'no zpe added'
	   elif ZPE_gas == 'vib' :
	       try: 	       
	           zpe=calc_zpe_vib(spec['vib'])
		   print 'vibrational zpe:',zpe
	       except KeyError: 
	           print >> check, "You don;t appeat to have any vibrations for species "+species.get(i)
		   print >> check, "zpe set to 0.0"
	           print >> check, "I would try the zero entropy and zpe surface model"
		   zpe=0	   
	   elif ZPE_gas == 'zpe' :
	       
	       try:      	       
	           zpe=spec['ZPE']
		   print 'pre_calculated zpe:',zpe
	       except KeyError:
	           print >> check, "you don't appear to have a ZPE for species "+species.get(i)
		   print >> check, "zpe set to 0.0"
	           print >> check, "I would try a different approximation or do some more DFT calcualtions"
		   zpe=0
	   		   		   	       	   
	          
           species_ZPE[i]=zpe
           
    return species_ZPE	 
	        
def get_p_w_species(species,datasurf=None,datagas=None,Temp=298., P_work=None):
	 #at the moment specific choices for each thermodynamic contribution need to be made (otherwise they are ignored and E is returned
    #At a later date coudl introduce default approximations too e.g. approx1=(blardy blah), approx2=another set of values
    
    species_pw=np.zeros(len(species))
    
    for i in range(0,len(species),1):
                   
       phase=species.get(i)[-1] 
       print >> check, 'Running for SPECIES ** :' , species.get(i)
              	           
       if phase == 'g':
           
           spec = datagas[0].get(datagas[1]).get(species.get(i))
           print 'spec',spec          	  
	          
	  #P_work
	   if P_work == None:
	       p_w = 0.0
	       print 'pressure work is:', p_w
	   else: 
	       try:
	           p_w=calc_p_work(Temp,P_work[species.get(i)])
		   print 'kTln(p/pstan):', p_w
	       except KeyError:
	           print >> check, "you appear to missing something todo with the pressure"
		   print >> check, "Check what you have entered for " + species.get(i)		   
		   p_w=0.0
           
    return species_pw   

def get_app_pot_species(species,datasurf=None,datagas=None,Temp=298., over_pot=None):
	 print 'to be implemented'
        
def calc_zpe_vib(FREQ):
    
        l = len(FREQ)
        dum =0.0        

        for i in xrange(0,int(l)):

                a = (h * c * L / 2000) * FREQ[i] / (96.485 )

#               print h * c * L / 2000

                dum = dum + a
                dum = dum

        ZPE = dum   
    
        return ZPE    	   

def calc_u_surf_vib(T, FREQ):
        #print FREQ
        l = len(FREQ)
        U_sum=0.0
        for i in xrange(0,int(l)):

            thetav = hc_over_k * FREQ[i]
        
            u = thetav / T
            a = u/(exp(u)-1.0)
           
            U = 8.3145 * (a)
	    
            U_sum += U*T
	    
        U_sum=U_sum/(96.485*1000.0)
	#print 'U_sum' , U_sum
        return U_sum

def calc_s_surf_vib(T,FREQ):
        #print FREQ
        l = len(FREQ)

        S_sum=0.0
        for i in xrange(0,int(l)):

            thetav = hc_over_k * FREQ[i]
        
            u = thetav / T
            a = u/(exp(u)-1.0)
            b = log(1.0-exp(-u))

            S = 8.3145 * (a-b)

            #print 'ess',S, FREQ, T
            if S > 8.3145*3. :  
                S = 8.3145*3
                #print 'ess_limit',S, FREQ, T

            S_sum += S

        S_sum = S_sum/(96.485*1000.0)
        #print 'S_sum', S_sum      
        return S_sum
   

def calc_h_gas_vib():
    print 'get_h_gas_shom not implemented'
    print 'setting h_gas to zero, I would advise using shomate parameters' 
    h_gas=0.0
    return h_gas
    
def calc_h_gas_shom(T, params):
        
	if len(params) == 8 :
	    print 'using shomate for single range'
	    data=params
	else:
	    for i in range(len(params)-1,-1,-1):
	        T1 = 0.0
		print len(params),'TESTS', i
	        if T <  params[i][0] and T > T1 :
		    data = params[i][1]
		    print 'Using shomate paramters for T = ',params[i][0]  
		    print 'which are : ' , params[i][0]
		    T1 = params[i][0]
	 
        #print data,'data',float(data[8])		
        #print data
        A = data[0]
        B = data[1]
        C = data[2]
        D = data[3]
        E = data[4]
        F = data[5]      
        H = data[7]

        T=T/1000.0
        Enthalpy = A*T + B*T*T/2. + C*T**3/3. + D*T**4/4. - (E/T) + F - H
        
        Enthalpy = Enthalpy/96.485
        #print Enthalpy
        return Enthalpy   
    
def calc_s_gas_shom(T,params):
        
        if len(params) == 8 :
	    print 'using shomate for single range'
	    data=params
	else:
	    for i in range(len(params)-1,-1,-1): 
	        T1 = 0.0
		print len(params),'TESTS', i
	        if T <  params[i][0] and T > T1 :	                  
		    data = params[i][1]
		    T1 = params[i][0]
        A = data[0]
#       print 'A',A
        B = data[1]
#       print B
        C = data[2]
#       print C
        D = data[3]
#       print D
        E = data[4]
#       print E
        G = data[6]
#       print 'G',G

        T=T/1000.0

        Entropy = A*log(T) + B*T + (C*pow(T,2.)/2.) + (D*pow(T,3.0)/3.) - (E/(2.*pow(T,2.0))) + G
        #print Entropy
        Entropy =Entropy/(96.485*1000.0)
        #print Entropy
        return Entropy    

def calc_s_gas_trans(T,M):
	
	S_trans = 1.5*R*log(M) + 2.5*R*log(T) - 1.1649*R
    
        return S_trans/(96.485*1000.0)
    
def calc_p_work(T,P):    
    
    pw=kbeV*T*log(P)

    return pw

def get_thermprop_reaction(species, thermprop_species, reactions, stoichiometry, act_bar=False):
    
    #this has been named thermprop, in anticipation of implementing S, H, etc
    thermprop_reaction=np.zeros(len(reactions))
    act_bar_reaction  =np.zeros(len(reactions))
    THP_LHS = np.zeros(len(reactions))
    THP_RHS = np.zeros(len(reactions))
    TS      = np.zeros(len(reactions))    

    for i in range(0,len(reactions),1):
        ky='r'+str(i)

	lhs_spec =  reactions[ky].get('lhs')
	lhs_stoich = stoichiometry[ky].get('lhs')
        #print lhs_spec
        for j in range(0,len(lhs_spec),1):
            for key in species.keys():
                #print 'LHS',key,species[key],lhs_spec[j], thermprop_species[key] 
                if species[key] == lhs_spec[j]:   
                   #print 'yaylhs'  
		   print 'lhs_stoich[j] ',lhs_stoich[j], len(lhs_spec), j,lhs_stoich               
                   THP_LHS[i] += thermprop_species[key]*lhs_stoich[j]      
    
        rhs_spec =  reactions[ky].get('rhs')
	rhs_stoich = stoichiometry[ky].get('rhs')
        #print rhs_spec
        for j in range(0,len(rhs_spec),1):
            for key in species.keys():
                #print 'RHS',key,species[key],rhs_spec[j], thermprop_species[key]
                if species[key] == rhs_spec[j]:
                   #print 'yayrhs'
                   THP_RHS[i] += thermprop_species[key]*rhs_stoich[j] 

        if act_bar == True:
            
            ts_spec =  reactions[ky].get('TS')
	    ts_stoich = stoichiometry[ky].get('TS')
            print 'ts+spec',ts_spec
            if ts_spec == None: 
                TS[i] += THP_RHS[i] 		
            else:
                for j in range(0,len(ts_spec),1):     
                    for key in species.keys():
                     #print 'TS',key,species[key],ts_spec[j], thermprop_species[key]
                        if species[key] == ts_spec[j]:
                            print 'yayts',thermprop_species[key], ts_stoich[j]
                            TS[i] += thermprop_species[key]*ts_stoich[j]
			    #if TS[i] <= THP_RHS[i]: TS[i] += THP_RHS[i]
			    #if TS[i] <= 0.0: TS[i] += 0.0
			     
                
        
        #print THP_LHS[i], THP_RHS[i]

    thermprop_reaction = THP_RHS - THP_LHS 
    act_bar_reaction = TS - THP_LHS
    for i in range(0,len(reactions),1):
        if act_bar_reaction[i] <= thermprop_reaction[i] : act_bar_reaction[i] = thermprop_reaction[i] 
        if act_bar_reaction[i] <= 0.0 : act_bar_reaction[i] = 0.0
    if act_bar == True:    
        return thermprop_reaction, act_bar_reaction
    else:
        return thermprop_reaction    

def calc_K_eqm(Temp, G_reac_array):
    K_eqm=np.zeros(len(G_reac_array))
    for i in range(0,len(G_reac_array),1):
        K_eqm[i]=exp(-G_reac_array[i]/(kbeV*Temp))
    return K_eqm   

def calc_k_kin(Temp, G_act_array):
    k_kin=np.zeros(len(G_act_array)) 
    for i in range(0,len(G_act_array),1):
        if G_act_array[i] <= 0.0 : G_act_array[i] = 0.0
        k_kin[i]=(kb*Temp/h)*exp(-G_act_array[i]/(kbeV*Temp))	
    return k_kin
    
def calc_k_kin2(Temp, H_act_array, S_act_array):
    k_kin=np.zeros(len(H_act_array)) 
    for i in range(0,len(H_act_array),1):
        if H_act_array[i] <= 0.0 : H_act_array[i] = 0.0
        k_kin[i]=(kb*Temp/h)*exp(S_act_array[i]/kbeV)*exp(-H_act_array[i]/(kbeV*Temp))	
    return k_kin    
    
def generate_stoichiometric_matrix(species, reactions, stoichiometries):
    stoich_mat=np.zeros((len(species)+1,len(reactions)),'d')
   
    for i in range(0,len(reactions),1):   
        #for j in range(0,len(species)+1,1):
            RHS = reactions['r'+str(i)].get('rhs')
            LHS = reactions['r'+str(i)].get('lhs')
	    RHS_stoic = stoichiometries['r'+str(i)].get('rhs')
	    LHS_stoic = stoichiometries['r'+str(i)].get('lhs')

            #print 'RHS_stoic' ,RHS_stoic
            #print 'LHS_stoic' , LHS_stoic

            #determine the number of empty sites and balance the sotichiometric *
            star_rhs=0
            star_lhs=0
            for j in range(0,len(species),1):
                for lenlhs in range(0,len(LHS),1):
                    if species[j] == LHS[lenlhs]:
                        #print 'species_lhs',j, species[j], LHS[lenlhs]
                        stoich_mat[j,i] += -LHS_stoic[lenlhs]
                        if species[j][-1] == '*':
                            star_lhs += LHS_stoic[lenlhs]
                for lenrhs in range(0,len(RHS),1):
                    if species[j] == RHS[lenrhs]:
                        #print 'species_rhs',j, species[j], RHS[lenrhs]
                        stoich_mat[j,i] += RHS_stoic[lenrhs]
                        if species[j][-1] == '*':
                            star_rhs += RHS_stoic[lenrhs]
            stoich_mat[-1,i] = star_lhs - star_rhs
                    
    #print 'STO', stoich_mat
    return stoich_mat

def calc_diff_eqns(cov, pres, stoich_mat, k_kin, species, remove_gas_species=True):

        rates, prespart = calc_elementary_rates(cov, pres, stoich_mat, k_kin)
        int = np.sum(rates*stoich_mat, axis = 1)
        short_species_list = species.copy()
        if remove_gas_species == True:
            for i in range(len(species)-1,-1,-1):
                
                if species[i][-1] == 'g' :
                        #np.delete(int,i,0) 
                        #del short_species_list[i] 
			int[i]=0        
                        #print ' deleted',species[i] 
	
        #int = np.sum(rates*stoich_mat, axis = 1)
        #print 'int', int
        if remove_gas_species == True:
            return int, short_species_list
        else:
            return int, pres
	    
def calc_diff_eqns2(pres_cov, stoich_mat, k_kin, species, fix_gas_species=False, flow=False, P_in= False):

        rates, prespart = calc_elementary_rates2(pres_cov, stoich_mat, k_kin)
        int = np.sum(rates*stoich_mat, axis = 1)
        short_species_list = species.copy()
        if fix_gas_species <> False:
           # for i in range(len(species)-1,-1,-1):
           #
	   for i in fix_gas_species:    
                if species[i][-1] == 'g' :
                        #np.delete(int,i,0) 
                        #del short_species_list[i] 
			int[i]=0        
                        #print ' deleted',species[i] 
	
	if flow <> False:
	    #print species
	    for spec in species:
	        if species[spec][-1] == 'g' :
	            int[spec] += (P_in[species[spec]]*flow)-(pres_cov[spec]*flow)
		    #print int[spec]
		    #print (P_in[species[spec]]*flow)-(pres_cov[spec]*flow)
		    #print 'check flow', P_in[species[spec]]*flow, pres_cov[spec]*flow
	            #print 'check flow',spec, P_in[species[spec]], spec, pres_cov[spec], flow 
	    #print '******************'
        #int = np.sum(rates*stoich_mat, axis = 1)
        #print 'int', int
        #if fix_gas_species == True:
        #    return int, pres_cov
        #else:
        #    return int, pres_cov	
	return int, pres_cov
	    
def calc_elementary_rates2(pres_cov, stoich_mat, k_kin):
    #calcualte the rates       
    #extract the lhs of the reaction 
    stoich_lhs=stoich_mat*(stoich_mat<0)
 
    #stoich_lhs=stoich_mat
    #print 'stoich_lhs', stoich_lhs
    #pres_cov[pres_cov < 0.0] = 1e-16
    #print ' pres_cov_calc',pres_cov

    prescovcomp = np.power(pres_cov,abs(np.transpose(stoich_lhs)))
    #print 'prescovcomp', prescovcomp 
   
    #get product of coverages (raised by their stoichiometric number)
    prescovpart = np.product(prescovcomp, axis=1)
    #print 'prescovpart',prescovpart        
    rates   = prescovpart*k_kin              
    #print 'RATES', rates

    return rates, prescovpart	        
	    
def calc_elementary_rates(cov, pres, stoich_mat, k_kin):
    #calcualte the rates       
    #extract the lhs of the reaction 
    stoich_lhs=stoich_mat*(stoich_mat<0)
    
    #raise the coverages and pressures to the power of the stoichiometric number
    #print 'cov',cov
    #print 'stoich_lhs',stoich_lhs
    #print np.transpose(stoich_lhs)
    covcomp = np.power(cov,abs(np.transpose(stoich_lhs)))
    
    #print 'covcomp',covcomp
    #print 'pres again', pres
    prescomp = np.power(pres,abs(np.transpose(stoich_lhs)))
    
    #print 'prescomp',prescomp
    
    #get product of coverages (raised by their stoichiometric number)
    covpart = np.product(covcomp, axis=1)
    
    #print 'covpart',covpart
    prespart= np.product(prescomp, axis=1)
    
    #print 'prespart',prespart
    #print 'k_kin',k_kin
    
    rates   = prespart*covpart*k_kin              

    return rates, prespart

def gen_pres_list(species, pressure_list):
    pres     = np.ones(len(species)+1)
    #print pres
    #print species
    for i in range(0,len(species),1):
        #print pressure_list.get(species[i])
        if pressure_list.get(species[i]) == None:
            pres[i]=1.0
        else:
            pres[i] = pressure_list.get(species[i])
            #print pres[i],pressure_list.get(species[i]) 
       #except KeyError:
       #     pres[i]=1.0
    pres[-1]=1.0
    #print pres
    return pres

def gen_cov_list(species, pressure_list, cov_init=None):
    if cov_init == None : cov_init=np.zeros(len(species))
    #if cov_init == None : cov_init= [1e-20]*len(species)
    cov     = np.ones(len(species)+1)
    #print cov
    cov_sum=0.0
    for i in range(0,len(species),1):
        if pressure_list.get(species[i]) == None:
            cov[i]=cov_init[i]
	    cov_sum += cov_init[i] 
        else:
            cov[i] = 1.0
            print cov[i],pressure_list.get(species[i]) 
       #except KeyError:
       #     pres[i]=1.0
    cov[-1]=1.0-cov_init[i]
    return cov

def create_reac_plot_data1(G_reaction,G_activation):

    from scipy.interpolate import InterpolatedUnivariateSpline, interp1d
    from scipy import interpolate    
    from pylab import plot
    step_space=1.
    step_width=0.5
    x=[]
    y=[]
    G=0.
    G_final=np.sum(G_reaction)

    nsteps = len(G_reaction)    
    
    x.append((0*step_space)-(step_width/2.))        
    x.append((0*step_space)+(step_width/2.))
    y.append(G) ; y.append(G)

    for i in range(1,nsteps+1,1):
        if G_activation[i-1] <> G_reaction[i-1]:
            #x.append(((i-1)*step_space)+(step_space/2.))
            #y.append(G+G_activation[i-1])
            fac=step_space-step_width

            #lhs_spline
            
            G_eff=-(G_reaction[i-1]-G_activation[i-1])
            x_lor = [x[-1]+0.1*fac,x[-1]+0.2*fac,x[-1]+0.3*fac,x[-1]+0.35*fac,x[-1]+0.5*fac,x[-1]+0.65*fac,x[-1]+0.7*fac,x[-1]+0.8*fac,x[-1]+0.9*fac]
            y_lor = [G+G_activation[i-1]*0.4, G+G_activation[i-1]*0.8 ,G+G_activation[i-1]*0.95, G+G_activation[i-1]*0.975, G+G_activation[i-1] \
                     , G+G_activation[i-1]-G_eff*0.025, G+G_activation[i-1]-G_eff*0.05, G+G_activation[i-1]-G_eff*0.2, G+G_activation[i-1]-G_eff*0.6]
            print x_lor, y_lor
            #plot(x_lor,y_lor)
            print str(987),x_lor
            spl=InterpolatedUnivariateSpline(x_lor,y_lor, k=3)
            x_ts=np.linspace(x[-1]+(fac*0.1), x[-1]+(fac*0.9), 21)
            y_ts=spl(x_ts)
            print str(12345),x_ts
            print y_ts
            for j in range(0,len(x_ts),1):
                x.append(x_ts[j]) ; y.append(y_ts[j])

            #x.append(((i-1)*step_space)+(step_space/2.)) ; y.append(G+G_activation[i-1])

        G=G+G_reaction[i-1]

        x.append((i*step_space)-(step_width/2.))        
        x.append((i*step_space)+(step_width/2.))
        y.append(G) ; y.append(G)  
        
    #plot(x,y)
    return x,y
    
def generate_labels_reac_plot(reac, stoich):
    
    labels=[]
    
    for i in range(0,len(reac),1):
        LHS = reac['r'+str(i)].get('lhs')
	#print 'LHS',LHS
	LHS_stoich = stoich['r'+str(i)].get('lhs')
	#print 'LHS_stoich',LHS_stoich
	lab=''  
	
        for j in range(0,len(LHS),1):
	    print j
	    print LHS_stoich[j]
	    if LHS_stoich[j] == 1 :
	        lab += LHS[j]+'+'
	    else:
	        lab += str(LHS_stoich[j])+LHS[j]+'+'	
	labels.append(lab[:-1])   
	     
    RHS = reac['r'+str(i)].get('rhs') 
    #print 'RHS', RHS
    RHS_stoich = stoich['r'+str(i)].get('rhs')      	
    #print RHS_stoich, RHS_stoich
    
    lab=''
    for j in range(0,len(RHS),1):
	if RHS_stoich[j] == 1 :
	    lab += RHS[j]+'+'
	else:
	    lab += str(RHS_stoich[j])+RHS[j]+'+'
    
    labels.append(lab[:-1])
    
    return labels	
	 
def generate_tex_labels_reac_plot(reac, stoich, tex=None):
    
    labels=[]
    
    from fractions import Fraction
            
    for i in range(0,len(reac),1):
        LHS = reac['r'+str(i)].get('lhs')
	#print 'LHS',LHS
	LHS_stoich = stoich['r'+str(i)].get('lhs')
	#print 'LHS_stoich',LHS_stoich
	  
	if tex is not None:
            lab='$\mathrm{'
        else:
            lab=''
        for j in range(0,len(LHS),1):
	    print j
	    print LHS_stoich[j]
	    #if LHS[j][-1] == 'g' and tex is not None :
	    
	    if tex is not None :
	       
	        if LHS_stoich[j] == 1 :	        
	            lab += tex[LHS[j]]+'+'
	        else:
	            lab += str(Fraction(LHS_stoich[j]))+tex[LHS[j]]+'+'
		    
            else:
	    
	        if LHS_stoich[j] == 1 :	        
	            lab += LHS[j]+'+'
	        else:
	            lab += str(Fraction(LHS_stoich[j]))+LHS[j]+'+'
	    		    		    		    
        if tex is not None:			
	    labels.append(lab[:-1]+'}$')  
	else:
	    labels.append(lab[:-1])
	print lab     
    RHS = reac['r'+str(i)].get('rhs') 
    #print 'RHS', RHS
    RHS_stoich = stoich['r'+str(i)].get('rhs')      	
    #print RHS_stoich, RHS_stoich
    if tex is not None:
        lab='$\mathrm{'
    else:
        lab=''
    for j in range(0,len(RHS),1):
        if tex is not None :
	    
	    if RHS_stoich[j] == 1 :	        
	        lab += tex[RHS[j]]+'+'
	    else:
	        lab += str(Fraction(RHS_stoich[j]))+tex[RHS[j]]+'+'
		    
        else:
	    
	    if RHS_stoich[j] == 1 :	        
	        lab += RHS[j]+'+'
	    else:
	        lab += str(Fraction(RHS_stoich[j]))+RHS[j]+'+'      
	
    if tex is not None:
        labels.append(lab[:-1]+'}$')
    else:
        labels.append(lab[:-1])
    print labels
    return labels

def create_2D_G_data_species(G_species, ref_species, species, alpha=True, TS=False, g_ref=None, labels=True):

    step_space=1.
    step_width=0.5
        
    ref_dic={}
    for sp in species:
        ref_dic[species.get(sp)]=G_species[len(ref_dic)]
    print ref_dic
    
    if g_ref == None:
        ref_values = calc_g_ref(G_species, ref_species, species)
    print 'ref_values ', ref_values
    
    x=[]
    y=[]
    sp_labels=[] 
    x_labels=[]   
    print 'starting'
    if alpha == True :
        for spec in sorted(species.values()):
	    if TS == False:	   
	        if spec[-3:-1] == 'TS' or spec[-1] == 'g' :
		    print 'skipping', spec  
		else:    
		    print spec, spec[-3:-1]  
                    parsed = parse(spec[:-1]).gen_elements()
		    g=0
		    dg=0
		    for i in parsed:
		        g += parsed[i]*ref_values[i]
			print 'Checking', g,parsed[i],ref_values[i]			 
	            print 'Checking', ref_dic[spec]		
		    dg = ref_dic[spec]-g
		    print 'Checking dg'	, dg
		    y.append(dg)	
		    y.append(dg)
		    sp_labels.append(spec)		    
                    print parsed
            if TS == True:
	        if spec[-1] == 'g' :
		    print 'skipping', spec  
		else:    
		    print spec, spec[-3:-1]  
                    parsed = parse(spec[:-1]).gen_elements()
		    g=0
		    dg=0
		    for i in parsed:
		        g += parsed[i]*ref_values[i]
			print g,parsed[i],ref_values[i]			 
	            print ref_dic[spec]		
		    dg = ref_dic[spec]-g
		    print dg	
		    y.append(dg)	
		    y.append(dg)
		    sp_labels.append(spec)		    
                    print parsed   		    
    
    for i in range(0,len(y)/2,1):
        x.append((i*step_space)-(step_width/2.))        
        x.append((i*step_space)+(step_width/2.))
	x_labels.append(i)              
        
    
    if alpha == False :
        print 'not implemented'
    
    if labels == True:	
        return x,y,sp_labels,x_labels
    else:
        return x,y		
    
def calc_g_ref(G_species, ref_species, species):
    
    #generate a list of species, dumstring can also be used to check the matrix for linalg is correct
    dumstring=''
    el_list={}
    for i in ref_species:
        dumstring = dumstring+i[:-2]
    el=parse(dumstring)	
    total_ref=el.gen_elements()
    #print 'total ref', total_ref	
    n_els=len(total_ref)    
    
    element_array=[]
    
    for i in ref_species:
    
        dum_array=np.zeros(n_els)
	#print 'IIII',i
        data = parse(i[:-2]).gen_elements()
	#print 'DATA',data
        for j in range(0,n_els,1):
	    for spec in data:
	        for n,d in enumerate(total_ref):
		    #print 'n,d',n,d
		    try:
	                if str(d) == spec:
		            dum_array[n] = data[spec] 
                    except KeyError:
		         print 'this species doesn;t contain element', str(d)
			 
	element_array.append(dum_array)

    #print element_array	
    #create g vector	
    G_vec=np.zeros(len(ref_species))
    #print G_species
    for i in range(0,len(G_vec),1):    
	for n in range(0,len(species),1):
	    try:
		if species[n] == ref_species[i]:
	            G_vec[i]=G_species[n]
	    except KeyError:
	        print 'This species is not there', ref_species[i]
    
    print 'G_vec',G_vec		
	
    #solve linear algebra - Im sure theres an easier way to tdo the above - 
    #I think I need to learn list comprehensions and regular expressions
    chem_pot = np.linalg.solve(element_array, G_vec)
    print 'Free energy reference values:'
    print 'Total reference:', total_ref
    print 'Chemical Potential:', chem_pot
    print '*****************************'
    print 'checking'
    print np.allclose(np.dot(element_array, chem_pot), G_vec)
    
    chem_pot_dic={}
    
    for n,d in enumerate(total_ref):
        chem_pot_dic[d]=chem_pot[len(chem_pot_dic)]
    
    #print chem_pot_dic
    
    return chem_pot_dic
    
def generate_tex_labels_species_plot(labels_uf, tex_labels):
    labels_f=[]

    for i in labels_uf:   
        labels_f.append('$'+tex_labels[i]+'$')
	
    return labels_f	
    
def cm1_to_eV(cm):
    eV = cm*c*L*h/(96.485*1000)
    return eV
      	   
	
def check_consistency(reac_dict, stoic_dict):
    
    print '####Running consistency check'
    print 'be warned, currently the parser only works with integers'
    print 'as such decimals will be rounded, so take care.'
    print 'Hacking the parser is on the todo list' 
    print '####'
    #jsut do for lhs and rhs for now, will look at ts later//////.....
   
    #generate a list of species, dumstring can also be used to check the matrix for linalg is correct
    for l in range(0,len(reac_dict),1):
        out = [[],[]]
	print '##Running for r'+str(l)
	print 'half reactions:'
        for step in ['lhs','rhs']:
            half_r=reac_dict['r'+str(l)][step]
            half_s=stoic_dict['r'+str(l)][step]	    	    	     
	    print step+' :',half_r, half_s
	    eldum=''	
	
            for m in range(0,len(half_r),1):
	        if half_r[m][-1] == 'g':
		    #print 'true'
	            #for n in range(0,len(half_s),1):  
	                eldum += '('+half_r[m][:-2]+')'+str(int(half_s[m]))       
	     
	        if half_r[m][-1] == '*':
		    #print 'also_true'
		    #for n in range(0,len(half_s),1):  
	                eldum += '('+half_r[m][:-1]+')'+str(int(half_s[m]))

            el = parse(eldum)
	    total_ref=el.gen_elements()
            print >> check, total_ref
    #    dumstring=''
    #    el_list={}
    #    for i in ref_species:
    #        dumstring = dumstring+i[:-2]
    #    el=parse(dumstring)	
    #    total_ref=el.gen_elements()
    #    #print 'total ref', total_ref	
    #    n_els=len(total_ref)    
    
    #element_array=[]
    exit
    
def gen_reac_dictionary(reactions):

    print "####Generating reaction and stoichiometric dictionaries"
    r_dic={}
    s_dic={}
    
    for i in range(0,len(reactions),1):
    
        lhs_dum=[]
	rhs_dum=[]
	ts_dum =[]
	lhs_sto=[]
	rhs_sto=[]
	ts_sto=[]
	
        spl=str.split(reactions[i])
	#print spl
	plus = [j for j,x in enumerate(spl) if x == '+']
	eq   = [j for j,x in enumerate(spl) if x == '=']
	act  = [j for j,x in enumerate(spl) if x == '\\TS']
	if len(act) == 0 : act.append(2000)
	print "plus,eq,act",plus, eq, act
	
	plus_lhs=[]
	plus_rhs=[]
	plus_ts=[]
	for k in plus:
	    if k < eq[0]:
	        plus_lhs.append(k)
	    if k > eq[0] and k < act[0]:	
	        plus_rhs.append(k)
	    if k > act[0]:
	        plus_ts.append(k)	
	#these loops do everythign to the left of any plus signs
	if len(plus_lhs) == 0:
	    lhs_sto.append(float(spl[0]))
	    lhs_dum.append(spl[1])
	elif len(plus_lhs) > 0 :
	    for k in plus_lhs:	        
	        lhs_dum.append(spl[k-1])
		lhs_sto.append(float(spl[k-2]))
	    lhs_dum.append(spl[plus_lhs[-1]+2])	
	    lhs_sto.append(float(spl[plus_lhs[-1]+1]))
	 	          
	if len(plus_rhs) == 0:
	    rhs_sto.append(float(spl[eq[0]+1]))
	    rhs_dum.append(spl[eq[0]+2])
	elif len(plus_rhs) > 0 :
	    for k in plus_rhs:	        
	        rhs_dum.append(spl[k-1])
		rhs_sto.append(float(spl[k-2]))
	    rhs_dum.append(spl[plus_rhs[-1]+2])	
	    rhs_sto.append(float(spl[plus_rhs[-1]+1]))
	    
	if len(plus_ts) == 0 and act[0] == 2000:
	    print "no transition state has been found r"+str(i)
	elif len(plus_ts) == 0 and act[0] < 2000: 
	    print "one TS species found", spl
	    ts_dum.append(spl[act[0]+2]) 
	    ts_sto.append(float(spl[act[0]+1])) 
	elif len(plus_ts) > 0 :
	    for k in plus_ts:	        
	        ts_dum.append(spl[k-1])
		ts_sto.append(float(spl[k-2]))
	    ts_dum.append(spl[plus_ts[-1]+2])	
	    ts_sto.append(float(spl[plus_ts[-1]+1]))
	    
        print '##stoichiometric dict for reac r'+str(i)	    
	print lhs_sto,rhs_sto, ts_sto
	print '##reaction dict for reac r'+str(i)	    
	print lhs_dum,rhs_dum, ts_dum 
				
	#now we have read in the reaction list, 
	#we need to put it into a dictinoary for further processing
	
	if act[0] == 2000:
	    r_dic['r'+str(i)]={'lhs':lhs_dum,'rhs':rhs_dum}
	    s_dic['r'+str(i)]={'lhs':lhs_sto,'rhs':rhs_sto}
	else:    
	    r_dic['r'+str(i)]={'lhs':lhs_dum,'rhs':rhs_dum,'TS':ts_dum}
	    s_dic['r'+str(i)]={'lhs':lhs_sto,'rhs':rhs_sto,'TS':ts_sto}
	    	
        print 'test', lhs_dum, rhs_dum, ts_dum, lhs_sto, rhs_sto, ts_sto
	
    return r_dic, s_dic  

	
    





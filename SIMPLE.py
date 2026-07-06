import numpy as np
import Testing
from scipy.linalg import solve_banded

def Simple(geo, var, Faces, Fields, Coupler, TDMA):

    Nx = geo.Nx
    Ny = geo.Ny


    #------------------ generating d section -----------------------------------------------------

    for i in range (Nx+1):
        for j in range(Ny):

            if Faces.a_P_u[i,j] < 1e-12:

                Coupler.d_we[i,j] = 0

            else:

                Coupler.d_we[i,j] = geo.deltay / Faces.a_P_u[i,j]

    for i in range (Nx):
        for j in range (Ny+1):

            if Faces.a_P_v[i,j] < 1e-12:

                 Coupler.d_ns[i,j] = 0

            else:

                Coupler.d_ns[i,j] = geo.deltax / Faces.a_P_v[i,j]

    
    #------------------- a' coefficient builder ----------------------------------------------------

    for i in range (Nx+1):
        for j in range (Ny):

            Coupler.a_EW_P[i,j] = var.density * Coupler.d_we[i,j] * geo.deltay


    for i in range (Nx):
        for j in range (Ny+1):

            Coupler.a_NS_P[i,j] = var.density * Coupler.d_ns[i,j] * geo.deltax

    for i in range(Nx):
        for j in range(Ny):

            Coupler.b[i,j] = var.density*(((Fields.u_psu[i,j]*geo.deltay) - (Fields.u_psu[i+1,j]*geo.deltay)) + ((Fields.v_psu[i,j]*geo.deltax) - (Fields.v_psu[i,j+1]*geo.deltax)))
    
    for i in range (Nx):
        for j in range (Ny):

            Coupler.a_P_P[i,j] = Coupler.a_EW_P[i,j] + Coupler.a_EW_P[i+1,j] + Coupler.a_NS_P[i,j] + Coupler.a_NS_P[i,j+1]

    
    #--------------------- TDMA P' --------------------------------
    # -------------------------------------------------------------

    NS_faces = Faces.NS_faces
    WE_faces = Faces.WE_faces
    upper = TDMA.upper_P
    diag = TDMA.diag_P
    lower = TDMA.lower_P
    P_prime = Fields.P_prime
    RHS = TDMA.RHS_P
   
    P_prime.fill(0)
   
    for i in range (Nx):

        #Resets the matrix when moving to a new i column
        upper.fill(0)
        diag.fill(0)
        lower.fill(0)
        RHS.fill(0)

        for j in range (Ny):
                
            # CASE 1 --------------------- fully embedded ------------------------------------
            
            if Faces.cells[i,j] == False:

                upper[j] = 0
                diag[j] = 1
                lower[j] = 0
                RHS[j] = 0


            # CASE 2---------------------- Checking for Adjacent Walls/Inlets  --------------------------------

            else:

                #Base case
                upper[j] = -Coupler.a_NS_P[i,j+1]
                diag[j] = Coupler.a_P_P[i,j]
                lower[j] = -Coupler.a_NS_P[i,j]
                RHS[j] = Coupler.b[i,j]

                # ---------------------- Checking North ------------------------


                if NS_faces[i,j+1] == 1:        # Wall,

                    upper[j] = 0
                    diag[j] -= Coupler.a_NS_P[i,j+1]

                elif NS_faces[i,j+1] == 2:      # Inlet,  zero-gradiant

                    upper[j] = 0
                    diag[j] -= Coupler.a_NS_P[i,j+1]

                elif NS_faces[i,j+1] == 3:       # Outlet, set value

                    upper[j] = 0
                    RHS[j] += Coupler.a_NS_P[i,j+1] * 0

                # ---------------------- Checking South ------------------------

                if NS_faces[i,j] == 1:

                    lower[j] = 0
                    diag[j] -= Coupler.a_NS_P[i,j]
                
                elif NS_faces[i,j] == 2:

                    lower[j] = 0
                    diag[j] -= Coupler.a_NS_P[i,j]


                elif NS_faces[i,j] == 3:

                    lower[j] = 0
                    RHS[j] += Coupler.a_NS_P[i,j] * 0
              

                # -------------------- Checking East ----------------------------

                if WE_faces[i+1,j] == 1:

                    diag[j] -= Coupler.a_EW_P[i+1,j]

                elif WE_faces[i+1,j] == 2:

                    diag[j] -= Coupler.a_EW_P[i+1,j]


                elif WE_faces[i+1,j] == 3:

                    RHS[j] += Coupler.a_EW_P[i+1,j] * 0

                else:                                           # 0 = Fluid

                    RHS[j] += Coupler.a_EW_P[i+1,j] * P_prime[i+1,j]


                # ------------------- Checking West --------------------------------

                if WE_faces[i,j] == 1:

                    diag[j] -= Coupler.a_EW_P[i,j]

                elif WE_faces[i,j] == 2:

                    diag[j] -= Coupler.a_EW_P[i,j]

                elif WE_faces[i,j] == 3:

                    RHS[j] += Coupler.a_EW_P[i,j] * 0

                else:                                           # 0 = Fluid 

                    if i > 0:

                        RHS[j] += Coupler.a_EW_P[i,j] * P_prime[i-1,j]

        
            """
            print('-------------------------------------------------------------------------')
            print('column', i)
            print('lower')
            print(np.array2string(lower, precision=2))
            print('diag')
            print(np.array2string(diag, precision=2))
            print('upper')
            print(np.array2string(upper, precision=2))
            print('RHS')
            print(np.array2string(RHS, precision=2))
            print('a_EW_P')
            print(np.array2string(np.flipud(Coupler.a_EW_P.T),formatter={'float_kind': lambda x: f"{x:10.3f}"}, max_line_width= 1000000000))
            print('a_NS_P')
            print(np.array2string(np.flipud(Coupler.a_NS_P.T),formatter={'float_kind': lambda x: f"{x:10.3f}"}, max_line_width= 1000000000))
            print('a_P_P')
            print(np.array2string(np.flipud(Coupler.a_P_P.T),formatter={'float_kind': lambda x: f"{x:10.3f}"}, max_line_width= 1000000000))
            """
       
        ab = np.zeros((3, geo.Ny))
        ab[0, 1:] = upper[:-1]
        ab[1,:] = diag
        ab[2,:-1] = lower[1:]
        Fields.P_prime[i,:] = solve_banded((1,1), ab, RHS)

    print('b')
    print(np.array2string(np.flipud(Coupler.b.T),formatter={'float_kind': lambda x: f"{x:10.3f}"}, max_line_width= 1000000000))
            

    #------------- Testing --------------------------
    Print = False
    if Print == True:
        Testing.P_prime(Fields)

        Print = False      
    #----------------------------------------------------



#=============================================================================================================
#================================ PRESSURE CORRECTION / RELAXATION ===========================================
#=============================================================================================================

    #Fields.u_old = np.copy(Fields.u)
    Fields.v_old = np.copy(Fields.v)

#========================================== Pressure =====================================================
    print(Fields.P)
    for i in range(geo.Nx):
        for j in range(geo.Ny):

            # ------- Relaxation ---------

            Fields.P[i,j] = Fields.P[i,j] + 0.3*Fields.P_prime[i,j] 

#======================================= U momentum =======================================================

    for i in range (Nx):
        for j in range(Ny-1):

            #Skipping Walls for pressure / relaxation 
            if Fields.u[i,j] == 0:

                continue

            else:

                # --------- Pressure correction ----------------

                Fields.u[i,j] = Fields.u_psu[i,j] + (Coupler.d_we[i,j]*(P_prime[i-1,j]-P_prime[i,j]))      

                # --------- Under-relaxation ------------------

                #West/East Inlet
                if Faces.WE_faces[i,j] == 2:                                            

                    #Fields.u[i,j] = var.u_inlet
                    continue

                #North Inlet
                elif Faces.NS_faces[i,j+1] == 2 or Faces.NS_faces[i-1,j+1] == 2:

                    continue

                #South Inlet
                elif Faces.NS_faces[i,j] == 2 or Faces.NS_faces[i-1,j] == 2:

                    continue
                
                # No inlet, normal relaxation correction
                else:

                    Fields.u[i,j] = 0.7*Fields.u[i,j] + (1 - 0.7)*Fields.u_old[i,j]   


    # ------------------ OUTLET COPYING -------------------------------------       
    # Algorithm is indentical to that of the one used in u_solver
     
    i = 0
    j = 0

    # West/East outlets

    for i in range (geo.Nx+1):
        for j in range (geo.Ny):

        
            if Faces.WE_faces[i,j] == 3:
        
                if var.Xdirec == 1:                                 # +x flow direc

                    Fields.u[i,j] = Fields.u[i-1,j]         # zero gradient pulled from west ( upstream )

                elif var.Xdirec == 0:                               # -x flow direx

                    Fields.u[i,j] = Fields.u[i+1,j]         # zero gradient pulled from east ( upstream )



    # North/South outlets 
    
    for i in range (geo.Nx):
        for j in range (geo.Ny):

            if Faces.WE_faces[i,j] != 1:                                      # 1 = Wall
                    
                #Checking for southern outlets
                if Faces.NS_faces[i,j] == 3 or Faces.NS_faces[i-1,j] == 3:

                    # -y flow direc
                    if var.Ydirec == 0:                                 # -y flow direc

                        Fields.u[i,j] = Fields.u[i,j+1]         # zero gradient pulled from north ( upstream )

                    # +y flow direc
                    elif var.Ydirec == 1:                                 

                        Fields.u[i,j] = Fields.u[i,j+1]         # zero gradient pulled from north 


                # Checking for nothern outlets
                if Faces.NS_faces[i,j+1] == 3 or Faces.NS_faces[i-1,j+1] == 3:
                    
                    # +y flow direc
                    if var.Ydirec == 1:                                

                        Fields.u[i,j] = Fields.u[i,j-1]         # zero gradient pulled from south ( upstream )

                    # -y flow direc
                    elif var.Ydirec == 0:                              

                        Fields.u[i,j] = Fields.u[i,j-1]         # zero gradient pulled from south ( upstream )



#================================= V Momentum =================================================================

    for i in range (Nx-1):
        for j in range(Ny):

            #Skipping Walls for pressure / relaxation
            if Fields.v[i,j] == 0:
                          
                continue

            else:

                # ----------- Pressure Correction ---------------------
                
                Fields.v[i,j] = Fields.v_psu[i,j] + (Coupler.d_ns[i,j]*(P_prime[i,j-1]-P_prime[i,j]))       # pressure correction

                #--------- Under-Relaxation ----------

                #North/South Inlet
                if Faces.NS_faces[i,j] == 2:

                    continue
                
                #East Inlet
                elif Faces.WE_faces[i,j] == 2 or Faces.WE_faces[i,j-1] == 2:

                    continue
                
                #West Inlet
                elif Faces.WE_faces[i+1,j] == 2 or Faces.WE_faces[i+1,j-1] == 2:

                    continue
                
                #No inlet, normal relaxation correction
                else:

                    Fields.v[i,j] = 0.7*Fields.v[i,j] + (1 - 0.7)*Fields.v_old[i,j]           


    # ------------------ OUTLET COPYING -------------------------------------       
    # Algorithm is indentical to that of the one used in v_solver

    i = 0
    j = 0

    # North/South outlets

    for i in range (geo.Nx):
        for j in range (geo.Ny+1):

        
            if Faces.NS_faces[i,j] == 3:
            
                
                if j == geo.Ny:
            
                    if var.Ydirec == 1:                                 # +y flow direc

                        Fields.v[i,j] = Fields.v[i,j-1]         # zero gradient pulled from south 

                    elif var.Ydirec == 0:                               # -y flow direx

                        Fields.v[i,j] = Fields.v[i,j-1]         # zero gradient pulled from south

                elif j == 0:

                    if var.Ydirec == 1:                                 # +y flow direc

                        Fields.v[i,j] = Fields.v[i,j+1]         # zero gradient pulled from south 

                    elif var.Ydirec == 0:                               # -y flow direx

                        Fields.v[i,j] = Fields.v[i,j+1]         # zero gradient pulled from south

                
    # West/East outlets 
    
    for i in range (geo.Nx):
        for j in range (geo.Ny):

            if Faces.NS_faces[i,j] != 1:                                      # 1 = Wall
                    
                #Checking for western outlets
                if Faces.WE_faces[i,j] == 3 or Faces.WE_faces[i,j-1] == 3:

                    # -x flow direc
                    if var.Xdirec == 0:                                   # -x flow direc

                        Fields.v[i,j] = Fields.v[i+1,j]         # zero gradient pulled from east ( upstream )

                    # +x flow direc
                    elif var.Xdirec == 1:                                 

                        Fields.v[i,j] = Fields.v[i+1,j]         # zero gradient pulled from east 


                # Checking for eastern outlets
                if Faces.WE_faces[i+1,j] == 3 or Faces.WE_faces[i+1,j-1] == 3:
                    
                    # +x flow direc
                    if var.Xdirec == 1:                                

                        Fields.v[i,j] = Fields.v[i-1,j]         # zero gradient pulled from west ( upstream )

                    # -x flow direc
                    elif var.Xdirec == 0:                              

                        Fields.v[i,j] = Fields.v[i-1,j]         # zero gradient pulled from west ( upstream )      


#===========================================================================================================================================            
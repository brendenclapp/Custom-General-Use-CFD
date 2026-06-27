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

            Coupler.b[i,j] = var.density*((Fields.u_psu[i,j]*geo.deltay) + (Fields.u_psu[i+1,j]*geo.deltay) + (Fields.v_psu[i,j]*geo.deltax) + (Fields.v_psu[i,j+1]*geo.deltax))
    
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
            
            if NS_faces[i,j] == 1 and NS_faces[i,j+1] == 1 and WE_faces[i,j] == 1 and WE_faces[i+1,j] == 1:

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

                    RHS[j] += Coupler.a_EW_P[i,j] * P_prime[i-1,j]

      
        
        ab = np.zeros((3, geo.Ny))
        ab[0, 1:] = upper[:-1]
        ab[1,:] = diag
        ab[2,:-1] = lower[1:]
        Fields.P_TDMA[i,:] = solve_banded((1,1), ab, RHS)


    print('lower')
    print(np.array2string(lower, precision=2))
    print('diag')
    print(np.array2string(diag, precision=2))
    print('upper')
    print(np.array2string(upper, precision=2))
    print('RHS')
    print(np.array2string(RHS, precision=2))

    #------------- Testing --------------------------
    Print = True
    if Print == True:
        Testing.P_TDMA(Fields)

        Print = False      
    #----------------------------------------------------




#======================================== CORRECTION ===================================================


    for i in range(geo.Nx):
        for j in range(geo.Ny):

            Fields.P[i,j] = Fields.P[i,j] + 0.5*Fields.P_prime[i,j] 


    for i in range (Nx):
        for j in range(Ny-1):

            Fields.u[i,j] = 0.5*Fields.u_psu[i,j] + 0.5*(Coupler.d_we[i,j]*(P_prime[i-1,j]-P_prime[i,j]))

    for i in range (Nx-1):
        for j in range(Ny):

            Fields.v[i,j] = 0.5*Fields.v_psu[i,j] + 0.5*(Coupler.d_ns[i,j]*(P_prime[i,j-1]-P_prime[i,j]))


    print('---------- final u--------------')
    print(np.array2string(np.flipud(Fields.u.T),formatter={'float_kind': lambda x: f"{x:6.3f}"}))
    print('---------- final v--------------')
    print(np.array2string(np.flipud(Fields.v.T),formatter={'float_kind': lambda x: f"{x:6.3f}"}))
    print('---------- final P--------------')
    print(np.array2string(np.flipud(Fields.P.T),formatter={'float_kind': lambda x: f"{x:6.3f}"}))
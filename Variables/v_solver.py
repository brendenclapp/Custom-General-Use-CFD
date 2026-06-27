
import numpy as np
import VanLeer
import Testing

def TVD_v(geo, var, Fields, Faces):

    i = 0
    j = 0
    Nx = geo.Nx
    Ny = geo.Ny


#----------------------------------------------------------------------------------------------------------------------------------------------------
#---------- Scalar centered face generation -------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

    # West/East Faces for the scalar CV
    for i in range (Nx+1):
        for j in range (Ny):                                                        # generates west/east faces, F = puA
            Faces.VF_we[i,j] = var.density * Fields.u[i,j] * geo.deltay

    # North/South Faces for the scalar CV
    for i in range (Nx):
        for j in range (Ny+1):                                                      # generates north/south faces, F = pvA
            Faces.VF_ns[i,j] = var.density * Fields.v[i,j] * geo.deltax


#----------------------------------------------------------------------------------------------------------------------------------------------------
# ------ Staggered V diffusion terms --------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

    # diffusion terms for the staggered grid are based off material / geometric properties so we can just calculate them straight up
   
    #Northern Southern Diffusion 
    for i in range (Nx):                                                              # in v staggered grid, v nodes will never be half distance from a NS wall                                              
        for j in range (Ny+1):

            if Fields.v[i,j] == 0:                                                      # when fluid == false, v[i,j] = 0, hence no diffusion term
                Faces.VD_ns[i,j] = 0

            else:
                Faces.VD_ns[i,j] = var.dyn_vis * geo.deltax / geo.deltay                # D = (μ/δ) * A

    
    #Western Diffusion 
    for i in range (1, Nx):                                                              # in a v staggered grid, v nodes will only ever be half distance (wall adjacent) from a WE wall 
        for j in range (Ny+1):

            if Fields.v[i,j] == 0:                                                      # if field is embedded into a wall
                Faces.VD_w[i,j] = 0                                                     # no diffusion


            # if its not embedded in a wall does the node have an embedded node to the west? 

            elif Fields.v[i-1,j] == 0:                                      # if it does have an embedded node to the west, then its either fully adjacent or half adjacent
   
                #Checking for western wall adjacencies

                if j == 0 or j == Ny:                                                   # if we are at the edges we just assume the block above extends past the border, so the edge never creates a half wall
                    
                    Faces.VD_w[i,j] = 2*(var.dyn_vis * geo.deltay / geo.deltax)           # D_fullwall = (μ/δ/2) * A  = 2D


                elif Fields.v[i-1,j-1] == 0 and Fields.v[i-1,j+1] == 0:                         # wall adjacent where entire western v face is wall adjacent
                    Faces.VD_w[i,j] = 2*(var.dyn_vis * geo.deltay / geo.deltax)                 # D_fullwall = (μ/δ/2) * A  = 2D

                else:                                                                           # wall adjacent where half of western v face is wall adjacent
                    Faces.VD_w[i,j] = 1.5*(var.dyn_vis * geo.deltay / geo.deltax)               # D_halfwall = ((μ/δ/2) * A/2) + (μ/δ * A/2) = 1.5D

            # if it doesnt have an embedded node directly west then it is fully interior

            else:                                                                               # northern face is entirely interior 
                Faces.VD_w[i,j] = var.dyn_vis * geo.deltay / geo.deltax                         # D_fullinterior = (μ/δ) * A


    #Eastern Diffusion 
    for i in range (Nx-1):                                                              # in a v staggered grid, v nodes will only ever be half distance (wall adjacent) from a WE wall 
        for j in range (Ny+1):

            # is the node embedded in a wall?
            if Fields.v[i,j] == 0:                                                      # if field is embedded into a wall
                Faces.VD_e[i,j] = 0                                                     # no diffusion


            # if its not embedded in a wall does the node have an embedded node to the east? 

            elif Fields.v[i+1,j] == 0:                                      # if it does have an embedded node to the east, then its either fully adjacent or half adjacent
   
                #Checking for western wall adjacencies

                if j == 0 or j == Ny:                                                           # if we are at the edges we just assume the block above extends past the border, so the edge never creates a half wall
                    Faces.VD_e[i,j] = 2*(var.dyn_vis * geo.deltay / geo.deltax)                 # D_fullwall = (μ/δ/2) * A  = 2D

                elif Fields.v[i+1,j-1] == 0 and Fields.v[i+1,j+1] == 0:                         # wall adjacent where entire eastern v face is wall adjacent
                    Faces.VD_e[i,j] = 2*(var.dyn_vis * geo.deltay / geo.deltax)                 # D_fullwall = (μ/δ/2) * A  = 2D

                else:                                                                           # wall adjacent where half of eastern v face is wall adjacent
                    Faces.VD_e[i,j] = 1.5*(var.dyn_vis * geo.deltay / geo.deltax)               # D_halfwall = ((μ/δ/2) * A/2) + (μ/δ * A/2) = 1.5D

            # if it doesnt have an embedded node directly west then it is fully interior

            else:                                                                               # eastern face is entirely interior 
                Faces.VD_e[i,j] = var.dyn_vis * geo.deltay / geo.deltax                         # D_fullinterior = (μ/δ) * A


    # Diffusion Face Field Print Statements
    """
    print('---------------VD_ns------------------')
    print(np.flipud(Faces.VD_ns.T))
    print('--------------VD_w--------------')
    print(np.flipud(Faces.VD_w.T))
    print('--------------VD_e--------------')
    print(np.flipud(Faces.VD_e.T))
    """

#----------------------------------------------------------------------------------------------------------------------------------------------------
#-------- staggered v faces, formed from averaging the scalar centered faces -------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

    for i in range (Nx):
        for j in range (1, Ny):
            
            # (1, Ny) no need to generate coefficients for j == 0 or j == Ny because the BC already tell us what we need to know
            # VF_we/ns already take into account the Mesh GUI via its reliance on v field

            # Staggered face values for v nodes formed from averaging the scalar centered faces

            Faces.Fn_v[i,j] = 0.5*(Faces.VF_ns[i,j]+Faces.VF_ns[i,j+1])
            Faces.Fs_v[i,j] = 0.5*(Faces.VF_ns[i,j]+Faces.VF_ns[i,j-1])
            Faces.Fe_v[i,j] = 0.5*(Faces.VF_we[i,j]+Faces.VF_we[i,j-1])
            Faces.Fw_v[i,j] = 0.5*(Faces.VF_we[i+1,j]+Faces.VF_we[i+1,j-1])

            # TVD neighbor coefficients formed from v CV faces and v CV diffusion terms

            Faces.a_w_v[i,j] = Faces.VD_w[i,j] + max(Faces.Fw_v[i,j], 0)
            Faces.a_e_v[i,j] = Faces.VD_w[i,j] + max(-Faces.Fe_v[i,j], 0)
            Faces.a_s_v[i,j] = Faces.VD_ns[i,j] + max(Faces.Fs_v[i,j], 0)
            Faces.a_n_v[i,j] = Faces.VD_ns[i,j] + max(-Faces.Fn_v[i,j],0)

            Faces.a_P_v[i,j] = Faces.a_w_v[i,j] + Faces.a_e_v[i,j] + Faces.a_s_v[i,j] + Faces.a_n_v[i,j] + (Faces.Fe_v[i,j]-Faces.Fw_v[i,j]) + (Faces.Fn_v[i,j]-Faces.Fs_v[i,j])

            Faces.dPdy[i,j] = (Fields.P[i,j]-Fields.P[i,j-1]) / geo.deltay                      # (Pp - Ps)/δ


            # TESTING ALGORITHM
            """
            print('a_w_v')
            print(np.flipud(Faces.a_w_v.T))
            """

            #INSERT TVD CORRECTION TERM ALGORITHM HERE
            VanLeer.Vanleer_v(i , j, var, geo, Fields, Faces)

    #print('SV_DC')
    #print(Faces.Sv_DC)


#===========================================================================================================|
# ----------------------------------------------TDMA -------------------------------------------------------|
#===========================================================================================================|

def TDMA_v(var, geo, Fields, Faces, TDMA):

    from scipy.linalg import solve_banded

    v = Fields.v
    NS_Faces = Faces.NS_faces
    upper = TDMA.upper_v
    diag = TDMA.diag_v
    lower = TDMA.lower_v
    RHS = TDMA.RHS_v

    for i in range (1, geo.Nx-1):

        #Resets the matrix when moving to a new i column
        upper.fill(0)
        diag.fill(0)
        lower.fill(0)
        RHS.fill(0)


        for j in range (1, geo.Ny):

            k = j-1


            # CASE 1 --------------------- fully embedded ------------------------------------
            if NS_Faces[i,j] == 1:

                upper[k] = 0
                diag[k] = 1
                lower[k] = 0
                RHS[k] = 0

            # CASE 2---------------------- Checking for Adjacent Walls/Inlets  --------------------------------
            elif NS_Faces[i,j] != 1:

                #Base case
                upper[k] = -Faces.a_n_v[i,j]
                diag[k] = Faces.a_P_v[i,j]
                lower[k] = -Faces.a_s_v[i,j]
                RHS[k] = (-Faces.dPdy[i,j]) + Faces.Sv_DC[i,j]       #adding source terms

                # ---------------------- Checking North ------------------------
                if NS_Faces[i,j+1] == 1:   # 1 = Wall

                    upper[k] = 0

                elif NS_Faces[i,j+1] == 2:  # 2 = inlet
                    
                    upper[k] = 0
                    RHS[k] += Faces.a_n_v[i,j]*v[i,j+1]
                
                elif NS_Faces[i,j+1] == 3:  # 3 = outlet

                    upper[k] = 0
                    RHS[k] += Faces.a_n_v[i,j]*v[i,j]               # vP = vN

                # ---------------------- Checking South ------------------------
  
                if NS_Faces[i,j-1] == 1:  # 1 = Wall

                    lower[k] = 0

                elif NS_Faces[i,j-1] == 2:  # 2 = inlet

                    lower[k] = 0
                    RHS[k] += Faces.a_s_v[i,j]*v[i,j-1]
                
                elif NS_Faces[i,j-1] == 3:    # 3 = outlet

                    lower[k] = 0
                    RHS[k] += Faces.a_s_v[i,j]*v[i,j]           # vP = vS

                # ---------------------- Checking West ---------------------------------------------

                if NS_Faces[i-1,j] == 1:                       # 1 = Wall

                    RHS[k] += 0

                elif NS_Faces[i-1,j] == 2:                     # 2 = inlet

                    RHS[k] += Faces.a_w_v[i,j]*v[i-1,j]

                elif NS_Faces[i-1,j] == 3:                     # 3 = outlet

                    RHS[k] += Faces.a_w_v[i,j]*v[i,j]
                
                else:                                           # 0 = Fluid

                    RHS[k] += Faces.a_w_v[i,j]*v[i-1,j]

                # ------------------------- Checking East ------------------------------------------

                if NS_Faces[i,j] == 1:

                    RHS[k] += 0

                elif NS_Faces[i,j] == 2:

                    RHS[k] += Faces.a_e_v[i,j]*v[i+1,j]

                elif NS_Faces[i,j] == 3:

                    RHS[k] += Faces.a_w_v[i,j]*v[i,j]

                else:                                           # 0 = Fluid

                    RHS[k] += Faces.a_e_v[i,j]*v[i+1,j]

        # MATRIX SOLVER ( in i loop, outside j loop)

        i_tdma = i -1

        ab = np.zeros((3, geo.Ny-1))
        ab[0, 1:] = upper[:-1]
        ab[1,:] = diag
        ab[2,:-1] = lower[1:]
        Fields.v_TDMA[i_tdma,:] = solve_banded((1,1), ab, RHS)

    #------------- Testing --------------------------
    Print = False
    if Print == True:
        Testing.v_TDMA(Fields)

        Print = False     
    #-----------------------------------------------


    # BOUNDARY ADDER

    for i in range (1, geo.Nx-1):
        for j in range(1, geo.Ny):

            Fields.v_psu[i,j] = Fields.v_TDMA[i-1,j-1]
    


    Fields.v_psu[:,0] = Fields.v[:,0]
    Fields.v_psu[:, geo.Ny] = Fields.v[:, geo.Ny]
    Fields.v_psu[0,:] = Fields.v[0,:]
    Fields.v_psu[geo.Nx-1,:] = Fields.v[geo.Nx-1,:]


    #------------------ Testing ------------------------
    Print = False
    if Print == True:
        Testing.v_psu_boundaries(Fields)

        Print = False                               


# ==================== OUTLET COPYING =========================

    i = 0
    j = 0

    # North/South outlets

    for i in range (geo.Nx-1):
        for j in range (geo.Ny+1):

        
            if NS_Faces[i,j] == 3:
            
                
                if j == geo.Ny:
            
                    if var.Ydirec == 1:                                 # +y flow direc

                        Fields.v_psu[i,j] = Fields.v_psu[i,j-1]         # zero gradient pulled from south 

                    elif var.Ydirec == 0:                               # -y flow direx

                        Fields.v_psu[i,j] = Fields.v_psu[i,j-1]         # zero gradient pulled from south

                elif j == 0:

                    if var.Ydirec == 1:                                 # +y flow direc

                        Fields.v_psu[i,j] = Fields.v_psu[i,j+1]         # zero gradient pulled from south 

                    elif var.Ydirec == 0:                               # -y flow direx

                        Fields.v_psu[i,j] = Fields.v_psu[i,j+1]         # zero gradient pulled from south

                
    # West/East outlets 
    
    for i in range (geo.Nx):
        for j in range (geo.Ny):

            if NS_Faces[i,j] != 1:                                      # 1 = Wall
                    
                #Checking for western outlets
                if Faces.WE_faces[i,j] == 3 or Faces.WE_faces[i,j-1] == 3:

                    # -x flow direc
                    if var.Xdirec == 0:                                   # -x flow direc

                        Fields.v_psu[i,j] = Fields.v_psu[i+1,j]         # zero gradient pulled from east ( upstream )

                    # +x flow direc
                    elif var.Xdirec == 1:                                 

                        Fields.v_psu[i,j] = Fields.v_psu[i+1,j]         # zero gradient pulled from east 


                # Checking for eastern outlets
                if Faces.WE_faces[i+1,j] == 3 or Faces.WE_faces[i+1,j-1] == 3:
                    
                    # +x flow direc
                    if var.Xdirec == 1:                                

                        Fields.v_psu[i,j] = Fields.v_psu[i-1,j]         # zero gradient pulled from west ( upstream )

                    # -x flow direc
                    elif var.Xdirec == 0:                              

                        Fields.v_psu[i,j] = Fields.v_psu[i-1,j]         # zero gradient pulled from west ( upstream )


    

    Print = False
    if Print == True:

        Testing.v_psu_final(Fields)

        Print = False
    

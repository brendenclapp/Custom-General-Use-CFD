import numpy as np
import VanLeer
import Testing

def TVD_u(geo, var, Fields, Faces):

    i = 0
    j = 0
    Nx = geo.Nx
    Ny = geo.Ny


#----------------------------------------------------------------------------------------------------------------------------------------------------
#---------- Scalar centered face generation -------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

    for i in range (Nx+1):                                                                  # generates west/east faces, F = puA
        for j in range (Ny):    
            Faces.UF_we[i,j] = var.density * Fields.u[i,j] * geo.deltay
     
    for i in range (Nx):                                                                    # generates north/south faces, F = pvA
        for j in range (Ny+1):
            Faces.UF_ns[i,j] = var.density * Fields.v[i,j] * geo.deltax


#----------------------------------------------------------------------------------------------------------------------------------------------------
# ------ Staggered U diffusion terms --------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

    # diffusion terms for the staggered grid are based off material / geometric properties so we can just calculate them straight up
   
    # Western Eastern Diffusion
    for i in range (Nx+1):                                # in u staggered grid, u nodes will never be half distance from a WE wall                                              
        for j in range (Ny):

            if Fields.u[i,j] == 0:                                                      # when fluid == false, u[i,j] = 0, hence no diffusion term
                Faces.UD_we[i,j] = 0

            else:
                Faces.UD_we[i,j] = ( var.dyn_vis * geo.deltay )/geo.deltax                # D = (μ/δ) * A

    #Northern Diffusion
    for i in range (Nx+1):                                  # in a u staggered grid, u nodes will only ever be half distance (wall adjacent) from a NS wall 
        for j in range (Ny-1):


            # is the node embedded in a wall?
            if Fields.u[i,j] == 0:                        # if node is fully embedded in wall
                Faces.UD_n[i,j] = 0                       # that node will have no northern diffusion


            # if its not embedded in a wall does the node have an embedded node above it? 

            elif Fields.u[i,j+1] == 0:                      # if it does have an embedded node above it, then its either fully adjacent or half adjacent

                #Checking for northern wall adjacencies

                if i == 0 or i == Nx:                                                   # if we are at the edges we just assume the block above extends past the border, so the edge never creates a half wall
                    Faces.UD_n[i,j] =  2*(var.dyn_vis * geo.deltax / geo.deltay)        # D_fullwall = (μ/δ/2) * A  = 2D

                elif Fields.u[i-1,j+1] == 0  and Fields.u[i+1,j+1] == 0:                # wall adjacent where entire northern u face is wall adjacent
                    Faces.UD_n[i,j] =  2*(var.dyn_vis * geo.deltax / geo.deltay)        # D_fullwall = (μ/δ/2) * A  = 2D

                else:                                                                   # wall adjacent where half of northern u face is wall adjacent
                    Faces.UD_n[i,j] = 1.5*(var.dyn_vis * geo.deltax / geo.deltay)       # D_halfwall = ((μ/δ/2) * A/2) + (μ/δ * A/2) = 1.5D

            # if it doesnt have an embedded node directly above it then it is fully interior

            else:                                                                       # northern face is entirely interior 
                Faces.UD_n[i,j] = var.dyn_vis * geo.deltax / geo.deltay                 # D_fullinterior = (μ/δ) * A

            

    #Southern Diffusion
    for i in range (Nx+1):                                  # in a u staggered grid, u nodes will only ever be half distance (wall adjacent) from a NS wall 
        for j in range (1, Ny):
                
            # is the node embedded in a wall?
            if Fields.u[i,j] == 0:                                                       # if node is fully embedded in wall
                Faces.UD_s[i,j] = 0                                                      # that node will have no southern diffusion


            # if its not embedded in a wall does the node have an embedded node below it? 

            elif  Fields.u[i,j-1] == 0:                                                 # if it does have an embedded node below it, then its either fully adjacent or half adjacent
                
                #Checking for southern wall adjacencies

                if i == 0 or i == Nx:                                                   # if we are at the edges we just assume the block above extends past the border, so the edge never creates a half wall
                    Faces.UD_s[i,j] =  2*(var.dyn_vis * geo.deltax / geo.deltay)        # D_fullwall = (μ/δ/2) * A  = 2D
                

                elif Fields.u[i-1,j-1] == 0 and Fields.u[i+1,j-1] == 0:                 # wall adjacent where entire southern u face is wall adjacent      
                    Faces.UD_s[i,j] =  2*(var.dyn_vis * geo.deltax / geo.deltay)        # D_fullwall = (μ/δ/2) * A  = 2D
                
                else:                                                                   # wall adjacent where half of southern u face is wall adjacent
                    Faces.UD_s[i,j] = 1.5*(var.dyn_vis * geo.deltax / geo.deltay)       # D_halfwall = ((μ/δ/2) * A/2) + (μ/δ * A/2) = 1.5D
                
            # if it doesnt have an embedded node directly below it then it is fully interior
            
            else:                                                                       # if southern face is entirely interior                                                                                          # D_fullwall = (μ/δ/2) * A  = 2D
                Faces.UD_s[i,j] = var.dyn_vis * geo.deltax / geo.deltay                  # D_fullinterior = (μ/δ) * A

    # Diffusion Face Field Print Statements
        
    print('--------------UD_we------------------')
    print(np.flipud(Faces.UD_we.T))
    print('--------------UD_n--------------')
    print(np.flipud(Faces.UD_n.T))
    print('--------------UD_s--------------')
    print(np.flipud(Faces.UD_s.T))
    

#----------------------------------------------------------------------------------------------------------------------------------------------------
#-------- staggered u faces, formed from averaging the scalar centered faces -------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

    for i in range (1, Nx):
        for j in range (Ny):
            
            # (1, Nx) no need to generate coefficients for i == 0 or i == Nx because the BC already tell us what we need to know
            # UF_we/ns already take into account the Mesh GUI via its reliance on u field

            # Staggered face values for u nodes formed from averaging the scalar centered faces
        
            Faces.Fw_u[i,j] = 0.5 * (Faces.UF_we[i,j] + Faces.UF_we[i-1,j])
            Faces.Fe_u[i,j] = 0.5 * (Faces.UF_we[i,j] + Faces.UF_we[i+1,j])
            Faces.Fs_u[i,j] = 0.5 * (Faces.UF_ns[i,j] + Faces.UF_ns[i-1,j])
            Faces.Fn_u[i,j] = 0.5 * (Faces.UF_ns[i,j+1] + Faces.UF_ns[i-1,j+1])

            # TVD neighbor coefficients formed from u CV faces and u CV diffusion terms
            Faces.a_w_u[i,j] = Faces.UD_we[i,j] + max(Faces.Fw_u[i,j], 0)
            Faces.a_e_u[i,j] = Faces.UD_we[i,j] + max(-Faces.Fe_u[i,j], 0)
            Faces.a_s_u[i,j] = Faces.UD_s[i,j] + max(Faces.Fs_u[i,j], 0)
            Faces.a_n_u[i,j] = Faces.UD_n[i,j] + max(-Faces.Fn_u[i,j],0)

            Faces.a_P_u[i,j] = Faces.a_w_u[i,j] + Faces.a_e_u[i,j] + Faces.a_s_u[i,j] + Faces.a_n_u[i,j] + (Faces.Fe_u[i,j]-Faces.Fw_u[i,j]) + (Faces.Fn_u[i,j]-Faces.Fs_u[i,j])

            Faces.dPdx[i,j] = (Fields.P[i,j]-Fields.P[i-1,j]) / geo.deltax                      # (Pp - Pw)/δ


            #INSERT TVD CORRECTION TERM ALGORITHM HERE
            VanLeer.Vanleer_u(i, j, var, geo, Fields, Faces)
    
    
    print('UD_we')
    print(np.array2string(np.flipud(Faces.UD_we.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('Fe_u')
    print(np.array2string(np.flipud(Faces.Fe_u.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('a_w_u')
    print(np.array2string(np.flipud(Faces.a_w_u.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('a_e_u')
    print(np.array2string(np.flipud(Faces.a_e_u.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('a_n_u')
    print(np.array2string(np.flipud(Faces.a_n_u.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('a_s_u')
    print(np.array2string(np.flipud(Faces.a_s_u.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('a_P_u')
    print(np.array2string(np.flipud(Faces.a_P_u.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    

#===========================================================================================================|
# ----------------------------------------------TDMA -------------------------------------------------------|
#===========================================================================================================|

def TDMA_u(var, geo, Fields, Faces, TDMA):

    from scipy.linalg import solve_banded

    u = Fields.u
    WE_Faces = Faces.WE_faces
    upper = TDMA.upper
    diag = TDMA.diag
    lower = TDMA.lower
    RHS = TDMA.RHS
    
    for i in range (1, geo.Nx):

        #Resets the matrix when moving to a new i column
        upper.fill(0)
        diag.fill(0)
        lower.fill(0)
        RHS.fill(0)

        for j in range (1, geo.Ny-1):

            k = j-1             # local index, so we can have range(1, geo.Ny-1) but still have the .TDMA's start at 0 index.
        
            # CASE 1 --------------------- fully embedded ------------------------------------

            if WE_Faces[i,j] == 1:         
                upper[k] = 0
                diag[k] = 1
                lower[k] = 0
                RHS[k] = 0


            # CASE 2---------------------- Checking for Adjacent Walls/Inlets  --------------------------------

            elif WE_Faces[i,j] != 1:
                
                #Base case
                upper[k] = -Faces.a_n_u[i,j]
                diag[k] = Faces.a_P_u[i,j]
                lower[k] = -Faces.a_s_u[i,j]
                RHS[k] = (-Faces.dPdx[i,j]) + Faces.Su_DC[i,j]       #adding source terms

                # ---------------------- Checking North ------------------------
                if WE_Faces[i,j+1] == 1:   # 1 = Wall

                    upper[k] = 0

                elif WE_Faces[i,j+1] == 2:  # 2 = inlet
                    
                    upper[k] = 0
                    RHS[k] += Faces.a_n_u[i,j]*u[i,j+1]
                
                elif WE_Faces[i,j+1] == 3:  # 3 = outlet

                    upper[k] = 0
                    RHS[k] += Faces.a_n_u[i,j]*u[i,j]               # uP = uN


                # ---------------------- Checking South ------------------------
  
                if WE_Faces[i,j-1] == 1:  # 1 = Wall

                    lower[k] = 0

                elif WE_Faces[i,j-1] == 2:  # 2 = inlet

                    lower[k] = 0
                    RHS[k] += Faces.a_s_u[i,j]*var.u_inlet
                
                elif WE_Faces[i,j-1] == 3:    # 3 = outlet

                    lower[k] = 0
                    RHS[k] += Faces.a_s_u[i,j]*u[i,j]           # uP = uS


                # ---------------------- Checking West ---------------------------------------------

                if WE_Faces[i,j] == 1:                          # 1 = Wall

                    RHS[k] += 0

                elif WE_Faces[i,j] == 2:                        # 2 = inlet

                    RHS[k] += Faces.a_w_u[i,j]*var.u_inlet

                elif WE_Faces[i,j] == 3:                        # 3 = outlet

                    RHS[k] += Faces.a_w_u[i,j]*u[i,j]
                
                else:                                           # 0 = Fluid

                    RHS[k] += Faces.a_w_u[i,j]*u[i-1,j]

                # ------------------------- Checking East ------------------------------------------

                if WE_Faces[i+1,j] == 1:                        # 1 = Wall

                    RHS[k] += 0

                elif WE_Faces[i+1,j] == 2:                      # 2 = inlet

                    RHS[k] += Faces.a_e_u[i,j]*var.u_inlet

                elif WE_Faces[i+1,j] == 3:                      # 3 = outlet

                    RHS[k] += Faces.a_e_u[i,j]*u[i,j]

                else:                                           # 0 = Fluid

                    RHS[k] += Faces.a_e_u[i,j]*u[i+1,j]
        

          
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
        
        # MATRIX SOLVER ( in i loop, outside j loop)

        i_tdma = i -1

        ab = np.zeros((3, geo.Ny-2))
        ab[0, 1:] = upper[:-1]
        ab[1,:] = diag
        ab[2,:-1] = lower[1:]
        Fields.u_TDMA[i_tdma,:] = solve_banded((1,1), ab, RHS)


    #------------- Testing --------------------------
    Print = True
    if Print == True:
        Testing.u_TDMA(Fields)

        Print = False      
    #----------------------------------------------------


    # BOUNDARY ADDER
    for i in range (1, geo.Nx):
        for j in range(1, geo.Ny-1):
            Fields.u_psu[i,j] = Fields.u_TDMA[i,j]


    Fields.u_psu[:,0] = Fields.u[:,0]
    Fields.u_psu[:, geo.Ny-1] = Fields.u[:, geo.Ny-1]
    Fields.u_psu[0,:] = Fields.u[0,:]
    Fields.u_psu[geo.Nx,:] = Fields.u[geo.Nx,:]

    #------------------ Testing ------------------------
    Print = True
    if Print == True:
        Testing.u_psu_boundaries(Fields)

        Print = False                               

#==================================== OUTLET COPYING =========================================================

    i = 0
    j = 0

    # West/East outlets

    for i in range (geo.Nx+1):
        for j in range (geo.Ny):

        
            if WE_Faces[i,j] == 3:
        
                if var.Xdirec == 1:                                 # +x flow direc

                    Fields.u_psu[i,j] = Fields.u_psu[i-1,j]         # zero gradient pulled from west ( upstream )

                elif var.Xdirec == 0:                               # -x flow direx

                    Fields.u_psu[i,j] = Fields.u_psu[i+1,j]         # zero gradient pulled from east ( upstream )



    # North/South outlets 
    
    for i in range (geo.Nx):
        for j in range (geo.Ny):

            if WE_Faces[i,j] != 1:                                      # 1 = Wall
                    
                #Checking for southern outlets
                if Faces.NS_faces[i,j] == 3 or Faces.NS_faces[i-1,j] == 3:

                    # -y flow direc
                    if var.Ydirec == 0:                                 # -y flow direc

                        Fields.u_psu[i,j] = Fields.u_psu[i,j+1]         # zero gradient pulled from north ( upstream )

                    # +y flow direc
                    elif var.Ydirec == 1:                                 

                        Fields.u_psu[i,j] = Fields.u_psu[i,j+1]         # zero gradient pulled from north 


                # Checking for nothern outlets
                if Faces.NS_faces[i,j+1] == 3 or Faces.NS_faces[i-1,j+1] == 3:
                    
                    # +y flow direc
                    if var.Ydirec == 1:                                

                        Fields.u_psu[i,j] = Fields.u_psu[i,j-1]         # zero gradient pulled from south ( upstream )

                    # -y flow direc
                    elif var.Ydirec == 0:                              

                        Fields.u_psu[i,j] = Fields.u_psu[i,j-1]         # zero gradient pulled from south ( upstream )

            


    # Print Statement ( 1 = PRINT, 0 = SKIP)

    Print = False
    if Print == True:

        Testing.u_psu_final(Fields)

        Print = False
    


def Vanleer_u(i , j, var, geo, Fields, Faces):
        
    # Makes use of the i, and j, from u_solver to solve for Su_DC in range of i:(1,Nx) and j(Ny)
    
    u = Fields.u

# X CASE 1 ------------------------------------        +x   flow       -------------------------------------------------------------------------------------

    if Faces.Fw_u[i,j] > 0  and Faces.Fe_u[i,j] > 0:

        var.Xdirec = 1
        

        # --------- r east ------------------

        
        if u[i,j] == 0 or u[i+1,j] == 0 or u[i-1,j] == 0:                            # if any of our stencil points are 0, then there is a wall there, and our stencil fails
                                                                                     # remember that u field adopted the geometric properties of the fluid TF field

            Faces.Su_DC[i,j] += 0                                                    # so our contribution from r_east is nothing

        else:                                                                        # if every stencil exists we can simply continue on with our calculations

            denom = (u[i+1,j]-u[i,j]) 
            if denom < 1e-12:                                                        # stops a division by zero if the denom of r is zero
               
                Faces.Su_DC[i,j] += 0

            else:                                                                    # if all stencil points valid, and no zero division, continue with VanLeer calculations

                r_east = (u[i,j]-u[i-1,j])/(u[i+1,j]-u[i,j])

                psi_east = (r_east + abs(r_east))/( 1 + r_east )

                Faces.Su_DC[i,j] += -0.5*Faces.Fe_u[i,j]*psi_east*(u[i+1,j]-u[i,j])

    
        # ---------- r west --------------------------------------------

        if u[i,j] == 0 or u[i-1,j] == 0 or u[i-2, j] == 0:                          # checking if all stencil points are valid

            Faces.Su_DC[i,j] += 0
        
        else:                                                                       # if all stencil points valid, and no zero division, continue with VanLeer calculations

            denom = (u[i,j]-u[i-1,j])
            if denom < 1e-12:                                                       # stops a division by zero if the denom of r is zero

                 Faces.Su_DC[i,j] += 0

            else:                                                                   # if all stencil points valid, and no zero division, continue with VanLeer calculations

                r_west = (u[i-1,j]-u[i-2,j])/(u[i,j]-u[i-1,j])

                psi_west = (r_west + abs(r_west))/(1 + r_west)

                Faces.Su_DC[i,j] += 0.5*Faces.Fw_u[i,j]*psi_west*(u[i,j]-u[i-1,j])



 # X CASE 2 ---------------------------------------------        -x  flow      ------------------------------------------------------------------

    elif Faces.Fw_u[i,j] < 0 and Faces.Fe_u[i.j] < 0:

        var.Xdirec = 0

        #------------ r east -----------------

        if u[i,j] == 0 or u[i+1, j] == 0 or u[i+2, j] == 0:                           # checking if all stencil points are valid

            Faces.Su_DC[i,j] += 0

        
        else:

            denom = (u[i+1,j]-u[i,j])
            if denom < 1e-12:                                                         # stops a division by zero if the denom of r is zero

                Faces.Su_DC[i,j] += 0

            else:                                                                     # if all stencil points valid, and no zero division, continue with VanLeer calculations

                r_east = (u[i+2,j]-u[i+1,j])/(u[i+1,j]-u[i,j])

                psi_east = (r_east + abs(r_east))/( 1 + r_east)
                
                Faces.Su_DC[i,j] += 0.5*Faces.Fe_u[i,j]*psi_east*(u[i+1,j]-u[i,j])



        # -----------  r west -------------------------

        if u[i,j] == 0 or u[i+1,j] == 0 or u[i-1,j] == 0:

            Faces.Su_DC[i,j] += 0

        else:


            denom = (u[i,j]-u[i-1,j])
            if denom < 1e-12:

                Faces.Su_DC[i,j] += 0

            else:

                r_west = (u[i+1,j]-u[i,j])/(u[i,j]-u[i-1,j])

                psi_west = (r_west + abs(r_west)) / (1 + r_west)

                Faces.Su_DC[i,j] += -0.5*Faces.Fw_u[i,j]*psi_west*(u[i,j]-u[i-1,j])


# Y CASE 1 ----------------------------------------------        +y       --------------------------------------------------------------------


    if Faces.Fn_u[i,j] > 0 and Faces.Fs_u[i,j] >0:

        var.Ydirec = 1


        # ------------------ r north -----------------------
        if j == geo.Ny-1:                                                             # bypass checking stencils if we are at the boundary edge
                                                                                    
            Faces.Su_DC[i,j] += 0

        elif u[i,j] == 0 or u[i,j-1] == 0 or u[i , j+1] == 0:                         # checking if all stencil points are valid

            Faces.Su_DC[i,j] += 0

        else:

            denom = (u[i,j+1]-u[i,j])
            if denom < 1e-12: 

                Faces.Su_DC[i,j] += 0

            else:

                r_north = (u[i,j]-u[i,j-1]) / (u[i,j+1]-u[i,j])

                psi_north = (r_north + abs(r_north))/(1 + r_north)

                Faces.Su_DC[i,j] += -0.5*Faces.Fn_u[i,j]*psi_north*(u[i,j+1]-u[i,j])

        

        # ------------------ r south -------------------------


        if u[i,j] == 0 or u[i,j-1] == 0 or u[i,j-2] == 0:                           # checking if all stencil points are valid

            Faces.Su_DC[i,j] += 0

        else:

            denom = (u[i,j]-u[i,j-1])
            if denom < 1e-12:

                Faces.Su_DC[i,j] += 0

            else:

                r_south = (u[i,j-1]-u[i,j-2])/(u[i,j]-u[i,j-1])

                psi_south = (r_south +abs(r_south))/( 1+ r_south)

                Faces.Su_DC[i,j] += 0.5*Faces.Fs_u[i,j]*psi_south*(u[i,j]-u[i,j-1]) 




# Y CASE 2 ------------------------------     -y       -------------------------------------------------------------

    elif Faces.Fn_u[i,j] < 0 and Faces.Fs_u[i,j] < 0:

        var.Ydirec = 0

        # ------------------ r north -----------------------

        if u[i,j] == 0 or u[i,j+1] == 0 or u[i,j+2] == 0:                     # checking if all stencil points are valid

            Faces.Su_DC[i,j] += 0

        
        else:

            denom = (u[i,j+1]-u[i,j])
            if denom <1e-12:
                    
                Faces.Su_DC[i,j] += 0

            else:

                r_north = (u[i,j+2] - u[i,j+1])/(u[i,j+1]-u[i,j])

                psi_north = (r_north + abs(r_north))/ ( 1 + r_north)

                Faces.Su_DC[i,j] += 0.5*Faces.Fn_u[i,j]*psi_north*(u[i,j+1]-u[i,j])


        # --------------------- r south ------------------------


        if u[i,j] == 0 or u[i,j+1] == 0 or u[i, j-1] == 0:

            Faces.Su_DC[i,j] += 0


        else:

            denom = ( u[i,j]-u[i,j-1])
            if denom < 1e-12:

                    Faces.Su_DC[i,j] += 0

            else:

                r_south = (u[i,j+1]-u[i,j])/( u[i,j]-u[i,j-1])

                psi_south = (r_south + abs(r_south))/( 1 + r_south)

                Faces.Su_DC[i,j] += -0.5*Faces.Fs_u[i,j]*psi_south*(u[i,j]-u[i,j-1])

#------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------


def Vanleer_v(i , j, var, geo, Fields, Faces):

    # Makes use of the i, and j, from u_solver to solve for Sv_DC 
    
    v = Fields.v

    # X CASE 1 ------------------------------------        +x   flow       -------------------------------------------------------------------------------------


    if Faces.Fw_v[i,j] > 0 and Faces.Fe_v[i,j] > 0:

        var.Xdirec = 1
        

        # --------- r east ------------------

        if i == geo.Nx-1:                                                           # bypass checking stencils if we are at the boundary edge
                                                                                    #kept just incase was important but dont think it is
            Faces.Sv_DC[i,j] += 0
        
        elif v[i,j] == 0 or v[i+1,j] == 0 or v[i-1,j] == 0:                            # if any of our stencil points are 0, then there is a wall there, and our stencil fails
                                                                                     # remember that u field adopted the geometric properties of the fluid TF field

            Faces.Sv_DC[i,j] += 0                                                    # so our contribution from r_east is nothing

        else:                                                                        # if every stencil exists we can simply continue on with our calculations

            denom = (v[i+1,j]-v[i,j]) 
            if denom < 1e-12:                                                        # stops a division by zero if the denom of r is zero
               
                Faces.Sv_DC[i,j] += 0

            else:                                                                    # if all stencil points valid, and no zero division, continue with VanLeer calculations

                r_east = (v[i,j]-v[i-1,j])/(v[i+1,j]-v[i,j])

                psi_east = (r_east + abs(r_east))/( 1 + r_east )

                Faces.Sv_DC[i,j] += -0.5*Faces.Fe_v[i,j]*psi_east*(v[i+1,j]-v[i,j])

        # ---------- r west --------------------------------------------

        if v[i,j] == 0 or v[i-1,j] == 0 or v[i-2, j] == 0:                          # checking if all stencil points are valid

            Faces.Sv_DC[i,j] += 0
        
        else:                                                                       # if all stencil points valid, and no zero division, continue with VanLeer calculations

            denom = (v[i,j]-v[i-1,j])
            if denom < 1e-12:                                                       # stops a division by zero if the denom of r is zero

                 Faces.Sv_DC[i,j] += 0

            else:                                                                   # if all stencil points valid, and no zero division, continue with VanLeer calculations

                r_west = (v[i-1,j]-v[i-2,j])/(v[i,j]-v[i-1,j])

                psi_west = (r_west + abs(r_west))/(1 + r_west)

                Faces.Sv_DC[i,j] += 0.5*Faces.Fw_v[i,j]*psi_west*(v[i,j]-v[i-1,j])



 # X CASE 2 ---------------------------------------------        -x  flow      ------------------------------------------------------------------

    elif Faces.Fw_v[i,j] < 0 and Faces.Fe_v[i,j] < 0:

        var.Xdirec = 0

        #------------ r east -----------------

        if v[i,j] == 0 or v[i+1, j] == 0 or v[i+2, j] == 0:                           # checking if all stencil points are valid

            Faces.Sv_DC[i,j] += 0

        
        else:

            denom = (v[i+1,j]-v[i,j])
            if denom < 1e-12:                                                         # stops a division by zero if the denom of r is zero

                Faces.Sv_DC[i,j] += 0

            else:                                                                     # if all stencil points valid, and no zero division, continue with VanLeer calculations

                r_east = (v[i+2,j]-v[i+1,j])/(v[i+1,j]-v[i,j])

                psi_east = (r_east + abs(r_east))/( 1 + r_east)
                
                Faces.Sv_DC[i,j] += 0.5*Faces.Fe_v[i,j]*psi_east*(v[i+1,j]-v[i,j])



        # -----------  r west -------------------------

        if v[i,j] == 0 or v[i+1,j] == 0 or v[i-1,j] == 0:

            Faces.Sv_DC[i,j] += 0

        else:


            denom = (v[i,j]-v[i-1,j])
            if denom < 1e-12:

                Faces.Sv_DC[i,j] += 0

            else:

                r_west = (v[i+1,j]-v[i,j])/(v[i,j]-v[i-1,j])

                psi_west = (r_west + abs(r_west)) / (1 + r_west)

                Faces.Sv_DC[i,j] += -0.5*Faces.Fw_v[i,j]*psi_west*(v[i,j]-v[i-1,j])


# Y CASE 1 ----------------------------------------------        +y       --------------------------------------------------------------------


    if Faces.Fn_v[i,j] > 0 and Faces.Fs_v[i,j] >0:

        var.Ydirec = 1


        # ------------------ r north -----------------------
        if v[i,j] == 0 or v[i,j-1] == 0 or v[i , j+1] == 0:                         # checking if all stencil points are valid

            Faces.Sv_DC[i,j] += 0

        else:

            denom = (v[i,j+1]-v[i,j])
            if denom < 1e-12: 

                Faces.Sv_DC[i,j] += 0

            else:

                r_north = (v[i,j]-v[i,j-1]) / (v[i,j+1]-v[i,j])

                psi_north = (r_north + abs(r_north))/(1 + r_north)

                Faces.Sv_DC[i,j] += -0.5*Faces.Fn_v[i,j]*psi_north*(v[i,j+1]-v[i,j])

        

        # ------------------ r south -------------------------


        if v[i,j] == 0 or v[i,j-1] == 0 or v[i,j-2] == 0:                           # checking if all stencil points are valid

            Faces.Sv_DC[i,j] += 0

        else:

            denom = (v[i,j]-v[i,j-1])
            if denom < 1e-12:

                Faces.Sv_DC[i,j] += 0

            else:

                r_south = (v[i,j-1]-v[i,j-2])/(v[i,j]-v[i,j-1])

                psi_south = (r_south +abs(r_south))/( 1+ r_south)

                Faces.Sv_DC[i,j] += 0.5*Faces.Fs_v[i,j]*psi_south*(v[i,j]-v[i,j-1]) 




# Y CASE 2 ------------------------------     -y       -------------------------------------------------------------

    elif Faces.Fn_v[i,j] < 0 and Faces.Fs_v[i,j] < 0:

        var.Ydirec = 0

        # ------------------ r north -----------------------

        if v[i,j] == 0 or v[i,j+1] == 0 or v[i,j+2] == 0:                     # checking if all stencil points are valid

            Faces.Sv_DC[i,j] += 0

        
        else:

            denom = (v[i,j+1]-v[i,j])
            if denom <1e-12:
                    
                Faces.Sv_DC[i,j] += 0

            else:

                r_north = (v[i,j+2] - v[i,j+1])/(v[i,j+1]-v[i,j])

                psi_north = (r_north + abs(r_north))/ ( 1 + r_north)

                Faces.Sv_DC[i,j] += 0.5*Faces.Fn_v[i,j]*psi_north*(v[i,j+1]-v[i,j])


        # --------------------- r south ------------------------


        if v[i,j] == 0 or v[i,j+1] == 0 or v[i, j-1] == 0:

            Faces.Sv_DC[i,j] += 0


        else:

            denom = ( v[i,j]-v[i,j-1])
            if denom < 1e-12:

                    Faces.Sv_DC[i,j] += 0

            else:

                r_south = (v[i,j+1]-v[i,j])/( v[i,j]-v[i,j-1])

                psi_south = (r_south + abs(r_south))/( 1 + r_south)

                Faces.Sv_DC[i,j] += -0.5*Faces.Fs_v[i,j]*psi_south*(v[i,j]-v[i,j-1])





def U_Field_adj(geo, var, Faces, Fields):

    for i in range(geo.Nx+1):
        for j in range(geo.Ny):

                if Faces.WE_faces[i,j] == 1:               # 1 = Wall
                     Fields.u[i,j] = 0

                elif Faces.WE_faces[i,j] == 2:             # 2 = Inlet
                     Fields.u[i,j] = var.u_inlet
                     


    return()

    
def V_Field_adj(geo, var, Faces, Fields):

    for i in range (geo.Nx):
        for j in range (geo.Ny+1):
             
            if Faces.NS_faces[i,j] == 1:                    # 1 = Wall
                  Fields.v[i,j] = 0
    

            elif Faces.NS_faces[i,j] == 2:                  # 2 = Inlet
                     Fields.v[i,j] = var.v_inlet

    

    return ()

     
    




def U_Field_adj(geo, var, Faces, Fields):

    for i in range(geo.Nx+1):
        for j in range(geo.Ny):

                if Faces.WE_faces[i,j] == 1:               # 1 = Wall
                    Fields.u[i,j] = 0
                
                elif Faces.WE_faces[i,j] == 2:             # 2 = WE Inlet
                    Fields.u[i,j] = var.u_inlet


    #North/South Inlets affecting u field
    for i in range(geo.Nx):
        for j in range(geo.Ny):
                
                if Faces.NS_faces[i,j+1] == 2 and Faces.NS_faces[i-1,j+1] == 2:
                    Fields.u[i,j] = var.u_inlet
                
                elif Faces.NS_faces[i,j] == 2 and Faces.NS_faces[i-1,j] == 2:
                     
                    Fields.u[i,j] = var.u_inlet
            
    


    return()

    
def V_Field_adj(geo, var, Faces, Fields):

    for i in range (geo.Nx):
        for j in range (geo.Ny+1):
             
            if Faces.NS_faces[i,j] == 1:                    # 1 = Wall
                  Fields.v[i,j] = 0
    

            elif Faces.NS_faces[i,j] == 2:                  # 2 = NS Inlet
                     Fields.v[i,j] = var.v_inlet


    #West/East Inlets affecting v field
    for i in range (geo.Nx):
        for j in range(geo.Ny):
             
            #West
            if Faces.WE_faces[i,j] == 2 and Faces.WE_faces[i,j-1] == 2:
                Fields.v[i,j] = var.v_inlet

            #East
            if Faces.WE_faces[i+1,j] == 2 and Faces.WE_faces[i+1, j-1] == 2:
                Fields.v[i,j] = var.v_inlet

    return ()

def P_field_adj(geo, Faces, Fields):
     
    for i in range(geo.Nx):
        for j in range(geo.Ny):
             
             if Faces.cells[i,j] == False:
                  
                  Fields.P[i,j] = 0


    return ()
    


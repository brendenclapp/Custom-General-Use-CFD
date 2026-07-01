import numpy as np

import Data
import TF_FieldAdj
import Variables.u_solver
import Variables.v_solver
import FACESGUI
import SIMPLE
import Testing


geo = Data.Geometry()
var = Data.Variables()
Fields = Data.Fields(geo)
Faces = Data.Faces(geo)
TDMA = Data.TDMA(geo)
Coupler = Data.Coupler(geo)


#-------- Runs GUI and generates cell type field -----------------

#TFgooey.CYCLE_GUI(Fields, geo)
FACESGUI.GUI(geo, Faces)

Print = False
if Print == True:
    Testing.geo(Faces)

    Print = False      

#=========================== LOOP START =================================================================

N = 1
i = 0
for i in range (N):
    #------- Adjust intial fields u/v in accordance to cell type field ---------------------------------------------------------------------------

    TF_FieldAdj.U_Field_adj(geo, var, Faces, Fields)                      # Adjusts the intial u-field in accordance with the cell_type geometry
    TF_FieldAdj.V_Field_adj(geo, var, Faces, Fields)                      # Adjust the intial v-field in accordance with the cell_type geometry
    TF_FieldAdj.P_field_adj(geo, Faces, Fields)

    

    print('----------u field----------------')
    print(np.flipud(Fields.u.T))
    print('----------v field----------------')
    print(np.flipud(Fields.v.T))
   
    #------------------- u variable solver ------------------------------------------

    Variables.u_solver.TVD_u(geo, var, Fields, Faces)
    Variables.u_solver.TDMA_u(var, geo, Fields, Faces, TDMA)
    


    #-------------------- v variable solver ------------------------------------------------

    Variables.v_solver.TVD_v(geo, var, Fields, Faces)
    Variables.v_solver.TDMA_v(var, geo, Fields, Faces, TDMA)

    #----------------------- COUPLER -------------------------------------------------------

    SIMPLE.Simple(geo, var, Faces, Fields, Coupler, TDMA)

    print('itt', i)
    #print('---------- final u--------------')
    #print(np.array2string(np.flipud(Fields.u.T),formatter={'float_kind': lambda x: f"{x:6.3f}"}))
    #print('---------- final v--------------')
    #print(np.array2string(np.flipud(Fields.v.T),formatter={'float_kind': lambda x: f"{x:6.3f}"}))
    print('---------- final P--------------')
    print(np.array2string(np.flipud(Fields.P.T),formatter={'float_kind': lambda x: f"{x:10.3f}"}, max_line_width= 1000000000))

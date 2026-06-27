import numpy as np

# the direct results of the TDMA algorithm, a version of u_psu that does not contain any boundary nodes
def u_TDMA(Fields):
    print('---------------TDMA u results, no boundary nodes---------------------------------')
    print(np.array2string(np.flipud(Fields.u_TDMA.T),formatter={'float_kind': lambda x: f"{x:6.3f}"}))
    print('\n')

# the direct results of the TDMA algorithm, a version of u_psu that does not contain any boundary nodes
def v_TDMA(Fields):
    print('---------------TDMA v results, no boundary nodes---------------------------------')
    print(np.array2string(np.flipud(Fields.v_TDMA.T),formatter={'float_kind': lambda x: f"{x:6.3f}"}))
    print('\n')

# version of u_psu where the boundary nodes have been added but outlets have not been updated
def u_psu_boundaries(Fields):
    print('---------------Boundary rows/columns added for u_psu ( outlets not updated yet ) ---------------------------------')
    print(np.array2string(np.flipud(Fields.u_psu.T),formatter={'float_kind': lambda x: f"{x:6.3f}"}))
    print('\n')

# version of v_psu where the boundary nodes have been added but outlets have not been updated
def v_psu_boundaries(Fields):
    print('---------------Boundary rows/columns added for v_psu ( outlets not updated yet ) ---------------------------------')
    print(np.array2string(np.flipud(Fields.v_psu.T),formatter={'float_kind': lambda x: f"{x:6.3f}"}))
    print('\n')

# Final variation of u_psu where all BCs have been implemented and zero gradiant outlets properly copied
def u_psu_final(Fields):
    
    print('---------------- final u_psu field -------------------------')
    print(np.array2string(np.flipud(Fields.u_psu.T),formatter={'float_kind': lambda x: f"{x:6.3f}"}))
    print('\n')


# Final variation of v_psu where all BCs have been implemented and zero gradiant outlets properly copied
def v_psu_final(Fields):
    
    print('---------------- final v_psu field -------------------------')
    print(np.array2string(np.flipud(Fields.v_psu.T),formatter={'float_kind': lambda x: f"{x:6.3f}"}))
    print('\n')


# Field that displays the P_prime correction field
def P_TDMA(Fields):

    print('---------------- P_prime field -------------------------')
    print(np.array2string(np.flipud(Fields.P_TDMA.T),formatter={'float_kind': lambda x: f"{x:6.3f}"}))
    print('\n')


def geo(Faces):

    print("Cells:")
    print("True = fluid, False = wall")
    print(np.flipud(Faces.cells.T))

    print("\nNS_faces / horizontal faces:")
    print("0 = fluid, 1 = wall, 2 = inlet, 3 = outlet")
    print(np.flipud(Faces.NS_faces.T))

    print("\nWE_faces / vertical faces:")
    print("0 = fluid, 1 = wall, 2 = inlet, 3 = outlet")
    print(np.flipud(Faces.WE_faces.T))
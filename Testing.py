import numpy as np

# the direct results of the TDMA algorithm, a version of u_psu that does not contain any boundary nodes
def u_TDMA(Fields):
    print('---------------TDMA u results, outlet not updated ---------------------------------')
    print(np.array2string(np.flipud(Fields.u.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('\n')

# the direct results of the TDMA algorithm, a version of u_psu that does not contain any boundary nodes
def v_TDMA(Fields):
    print('---------------TDMA v results, no boundary nodes---------------------------------')
    print(np.array2string(np.flipud(Fields.v_TDMA.T),formatter={'float_kind': lambda x: f"{x:6.3f}"}))
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
def P_prime(Fields):

    print('---------------- P_prime field -------------------------')
    print(np.array2string(np.flipud(Fields.P_prime.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
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

def a_u(Faces):

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
    

def TDMA_Testing(TDMA, i):

    print('-------------------------------------------------------------------------')
    print('column', i)
    print('lower')
    print(np.array2string(TDMA.lower, precision=6))
    print('diag')
    print(np.array2string(TDMA.diag, precision=6))
    print('upper')
    print(np.array2string(TDMA.upper, precision=6))
    print('RHS')
    print(np.array2string(TDMA.RHS, precision=6))
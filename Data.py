import numpy as np

class Geometry:

    def __init__(self):

        self.Nx = 8                                         # number of cells x direction
        self.Ny = 8                                         # number of cells in y direction
        self.deltax = 0.4                                   # length of unit cell (m)
        self.deltay = 0.1                                   # height of unit cell  (m)

class Variables:

    def __init__(self):
        self.density = 998                                  # density of fluid (1000kg/m^3)
        self.dyn_vis = 1e-6                                 # dynamic viscocity of fluid (kPa*s)
        self.u_inlet = 5                                    # u velocity of fluid coming into the system
        self.v_inlet = 0.5                                  # v velocity of fluid coming into the system
        self.Xdirec = 0                                     # boolean switch for fluid moving +x or -x
        self.Ydirec = 0                                     # boolean switch for fluid moving +y or -y
        self.P_anchor = 0                                   # Pressure anchor for pressure solver

class Fields:        

    def __init__(self, geo):
        self.u = np.full((geo.Nx+1, geo.Ny), 3.0)                # initial u_field with value set to 3 rn of u velocity
        self.u_old = np.full((geo.Nx+1, geo.Ny), 3)            # previous itteration of the u_field
        self.v = np.full((geo.Nx, geo.Ny+1), 0.3)              # initial v_field with value set to 1 rn of v velocity
        self.v_old = np.full((geo.Nx, geo.Ny+1), 0.3)          # previous itteration of the v_field
        self.P = np.full((geo.Nx, geo.Ny), 0.0)                  # Pressure field (kPa)
        self.P_prime = np.zeros((geo.Nx, geo.Ny))              # Pressure correction field
        self.u_psu = np.zeros((geo.Nx+1, geo.Ny))              # field that holds psuedo u, the TDMA u solution, which has not yet been pressure corrected
        self.u_TDMA = np.zeros((geo.Nx-1, geo.Ny-2))           # field that holds the smaller column matrix of TDMA results before placed into u_psu 
        self.v_psu = np.zeros((geo.Nx, geo.Ny+1))              # field that holds psuedo v, the TDMA v solution, which has not yet been pressure corrected
        self.v_TDMA = np.zeros((geo.Nx-2, geo.Ny-1))           # field that holds the smaller column matrix of TDMA results before placed into v_psu 

class Faces:

    def __init__(self, geo):

        # F reprents the convective mass flux per unit area at the cell face ( convection )
        # D represent the diffusion conductance at the cell face    ( diffusion )

        # ------------ U SOLVER ------------------------------------------------------------------------------------------------------------

        self.UF_we = np.zeros((geo.Nx+1, geo.Ny))           # West East faces for the scalar CV, specifically for the u variable / solver
        self.UF_ns = np.zeros((geo.Nx, geo.Ny+1))           # North South faces for the scalar CV, specifically for the u variable / solver

        self.UD_we = np.zeros((geo.Nx+1, geo.Ny))           # West East Diffusion values for the u CV, specifically for the u variable / solver
        self.UD_n = np.zeros((geo.Nx+1, geo.Ny))            # North Diffusion values for the u CV, specifically for the u variable / solver
        self.UD_s = np.zeros((geo.Nx+1, geo.Ny))            # South Diffusion values for the u CV, specifically for the u variable / solver

        self.Fw_u = np.zeros((geo.Nx+1, geo.Ny))            # Convection Flux for western u CV 
        self.Fe_u = np.zeros((geo.Nx+1, geo.Ny))            # Convection Flux for eastern u CV 
        self.Fn_u = np.zeros((geo.Nx+1, geo.Ny))            # Convection Flux for northern u CV 
        self.Fs_u = np.zeros((geo.Nx+1, geo.Ny))            # Convection Flux for southern u CV 

        self.a_w_u = np.zeros((geo.Nx+1, geo.Ny))           # western variable coefficient of u CV for U solver
        self.a_e_u = np.zeros_like(self.a_w_u)              # eastern variable coefficient of u CV for U solver
        self.a_n_u = np.zeros_like(self.a_w_u)              # northern variable coefficient of u CV for U solver
        self.a_s_u = np.zeros_like(self.a_w_u)              # southern variable coefficient of u CV for U solver
        self.a_P_u = np.zeros_like(self.a_w_u)              # center variable coefficient of u CV for U solver

        self.dPdx = np.zeros_like(self.a_w_u)               # horizontal pressure gradiant located at each ( shared w/ ) U node

        self.Su_DC = np.zeros((geo.Nx+1, geo.Ny))           # TVD Correction source term for u solver

        self.NS_faces = np.zeros((geo.Nx, geo.Ny + 1))      # horizontal faces for GUI
        self.WE_faces = np.zeros((geo.Nx + 1, geo.Ny))      # vertical faces for GUI
        self.cells = np.ones((geo.Nx, geo.Ny), dtype=bool)  # T/F for fluid or walls for scalar CV

        # ------------ V SOLVER ------------------------------------------------------------------------------------------------------------

        self.VF_we = np.zeros((geo.Nx+1, geo.Ny))           # West East faces for the scalar CV, specifically for the v variable / solver
        self.VF_ns = np.zeros((geo.Nx, geo.Ny+1))           # North South faces for the scalar CV, specifically for the v variable / solver

        self.VD_ns = np.zeros((geo.Nx, geo.Ny+1))           # North south Diffusion values for the v CV, specifically for the v variable / solver        
        self.VD_w = np.zeros((geo.Nx, geo.Ny+1))            # West Difussion values for the v CV, specifically for the v variable / solver
        self.VD_e = np.zeros((geo.Nx, geo.Ny+1))            # East Difussion values for the v CV, specifically for the v variable / solver

        self.Fw_v = np.zeros((geo.Nx, geo.Ny+1))            # Convection Flux for western v CV 
        self.Fe_v = np.zeros((geo.Nx, geo.Ny+1))            # Convection Flux for eastern v CV 
        self.Fn_v = np.zeros((geo.Nx, geo.Ny+1))            # Convection Flux for northern v CV 
        self.Fs_v = np.zeros((geo.Nx, geo.Ny+1))            # Convection Flux for southern v CV 

        self.a_w_v = np.zeros((geo.Nx, geo.Ny+1))           # western variable coefficient of v CV for V solver
        self.a_e_v = np.zeros_like((self.a_w_v))            # eastern variable coefficient of v CV for V solver
        self.a_n_v = np.zeros_like(self.a_w_v)              # northern variable coefficient of v CV for V solver
        self.a_s_v = np.zeros_like(self.a_w_v)              # southern variable coefficient of v CV for V solver
        self.a_P_v = np.zeros_like(self.a_w_v)              # center variable coefficient of v CV for V solver

        self.dPdy = np.zeros_like(self.a_w_v)               # vertical pressure gradiant located at each ( shared w/ ) V node

        self.Sv_DC = np.zeros((geo.Nx, geo.Ny+1))           # TVD Correction source term for v solver

class TDMA:

    def __init__(self, geo):

        self.upper = np.zeros((geo.Ny-2))
        self.diag = np.zeros((geo.Ny-2))
        self.lower = np.zeros((geo.Ny-2))
        self.RHS = np.zeros((geo.Ny-2))

        self.upper_v = np.zeros((geo.Ny-1))
        self.diag_v = np.zeros((geo.Ny-1))
        self.lower_v = np.zeros((geo.Ny-1))
        self.RHS_v = np.zeros((geo.Ny-1))

        self.upper_P = np.zeros((geo.Ny))
        self.diag_P = np.zeros((geo.Ny))
        self.lower_P = np.zeros((geo.Ny))
        self.RHS_P = np.zeros((geo.Ny))

class Coupler:

    def __init__(self,geo):

        self.d_we = np.zeros((geo.Nx+1, geo.Ny))
        self.d_ns = np.zeros((geo.Nx, geo.Ny+1))

        self.a_EW_P = np.zeros((geo.Nx+1, geo.Ny))
        self.a_NS_P = np.zeros((geo.Nx, geo.Ny+1))
        self.a_P_P = np.zeros((geo.Nx, geo.Ny))

        self.b = np.zeros((geo.Nx, geo.Ny))


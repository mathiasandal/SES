import numpy as np

'''Contains all functions related to things retrieved from VERES'''


def read_re7_file(filename):
    # TODO: Fix/add documentation
    """
    Reads a *.re7 file created by VERES

    A description of a *.re7 file is found in page 114 of the ShipX Vessel Responses (VERES) User's manual.
    """

    f = open(filename, "r")  # open file for reading

    for i in range(6):  # skips first six lines
        f.readline()

    # Read in parameters
    RHOSW, GRAV = [float(i) for i in f.readline().split()]
    LPP, BREADTH, DRAUGHT = [float(i) for i in f.readline().split()]
    LCG, VCG = [float(i) for i in f.readline().split()]

    # saves line looking like "-1 FILEVER (=2)" in Veres_Manual.pdf
    FILEVER = [int(i) for i in f.readline().split()]

    # Read in number of velocities, headings, frequencies and degrees of freedom
    NOVEL, NOHEAD, NOFREQ, NDOF = [int(i) for i in f.readline().split()]

    # Initialize vectors to contain velocities, headings and frequencies
    VEL = np.zeros(NOVEL)  # [m/s] Velocity of the vessel
    HEAD = np.zeros(NOHEAD)  # [deg] Wave heading
    FREQ = np.zeros(NOFREQ)  # [rad/s] Wave frequency

    # Initialize values to be stores for each velocity
    SINK = np.zeros(NOVEL)  # [m] Sink
    TRIM = np.zeros(NOVEL)  # [deg] Trim
    XMTN = np.zeros(NOVEL)  # [m] X-pos. of the motion coordinate system (relative to Lpp/2)
    ZMTN = np.zeros(NOVEL)  # [m] Z-pos. of the mo􀆟on coordinate system (relative to BL)

    # Read in mass matrix
    VMAS = np.zeros([NDOF, NDOF])  # Initialize mass matrix

    for j in range(NDOF):  # Read and converts each line from string to floats
        VMAS[j, :] = [float(i) for i in f.readline().split()]

    # Have one (NDOF x NDOF) matrix for each combinations of velocity, heading and wave frequency
    ADDMAS = np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Hydrodynamic added mass
    ADDADDMAS = np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Additional added mass
    DAMP = np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Hydrodynamic damping
    ADDDAMP = np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Additional damping
    REST = np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Restoring
    ADDREST = np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Additional restoring

    # Initialize viscous roll damping
    VISCDL = np.zeros([NOVEL, NOHEAD, NOFREQ])  # Viscous roll damping, linear part
    VISCDN = np.zeros([NOVEL, NOHEAD, NOFREQ])  # Viscous roll damping, nonlinear part
    VISCDNL = np.zeros([NOVEL, NOHEAD, NOFREQ])  # Viscous roll damping, linearized

    for i in range(NOVEL):
        VEL[i], SINK[i], TRIM[i], XMTN[i], ZMTN[i] = [float(m) for m in f.readline().split()]
        for j in range(NOHEAD):
            HEAD[j] = float(f.readline())  # Should only contain one element
            for k in range(NOFREQ):
                FREQ[k] = float(f.readline())  # Should only contain one element

                for m in range(NDOF):
                    ADDMAS[i, j, k, m, :] = [float(i) for i in f.readline().split()]

                for m in range(NDOF):
                    ADDADDMAS[i, j, k, m, :] = [float(i) for i in f.readline().split()]

                for m in range(NDOF):
                    DAMP[i, j, k, m, :] = [float(i) for i in f.readline().split()]

                for m in range(NDOF):
                    ADDDAMP[i, j, k, m, :] = [float(i) for i in f.readline().split()]

                for m in range(NDOF):
                    REST[i, j, k, m, :] = [float(i) for i in f.readline().split()]

                for m in range(NDOF):
                    ADDREST[i, j, k, m, :] = [float(i) for i in f.readline().split()]

                # Read in viscous roll damping
                VISCDL[i, j, k], VISCDN[i, j, k], VISCDNL[i, j, k] = [float(i) for i in f.readline().split()]
    return VMAS, ADDMAS, DAMP, REST, VEL, HEAD, FREQ, XMTN, ZMTN


def read_re8_file(filename):
    # TODO: Fix/add documentation
    """
    Reads a *.re8 file created by VERES

    A description of a *.re8 file is found in page 116 of the ShipX Vessel Responses (VERES) User's manual.

    Only works for High-Speed formulation in VERES
    """

    f = open(filename, 'r')

    for i in range(6):  # skips first six lines
        f.readline()

    # Read in parameters
    RHOSW, GRAV = [float(i) for i in f.readline().split()]
    LPP, BREADTH, DRAUGHT = [float(i) for i in f.readline().split()]
    LCG, VCG = [float(i) for i in f.readline().split()]

    # Read in number of velocities, headings, frequencies and degrees of freedom
    NOVEL, NOHEAD, NOFREQ, NDOF = [int(i) for i in f.readline().split()]


    # Initialize vectors to contain velocities, headings and frequencies
    VEL = np.zeros(NOVEL)  # [m/s] Velocity of the vessel
    HEAD = np.zeros(NOHEAD)  # [deg] Wave heading
    FREQ = np.zeros(NOFREQ)  # [rad/s] Wave frequency

    # Initialize values to be stores for each velocity
    SINK = np.zeros(NOVEL)  # [m] Sink
    TRIM = np.zeros(NOVEL)  # [deg] Trim
    XMTN = np.zeros(NOVEL)  # [m] X-pos. of the motion coordinate system (relative to Lpp/2)
    ZMTN = np.zeros(NOVEL)  # [m] Z-pos. of the mo􀆟on coordinate system (relative to BL)

    # Initialize Force components
    # Real parts
    REFORCE = np.zeros([NDOF, NOFREQ, NOHEAD, NOVEL])  # Real part
    IMFORCE = np.zeros([NDOF, NOFREQ, NOHEAD, NOVEL])  # Imaginary part

    for i in range(NOVEL):
        VEL[i], SINK[i], TRIM[i], XMTN[i], ZMTN[i] = [float(m) for m in f.readline().split()]
        for j in range(NOHEAD):
            HEAD[j] = float(f.readline())  # Should only contain one element
            for k in range(NOFREQ):
                FREQ[k] = float(f.readline())  # Should only contain one element

                for m in range(NDOF):
                    REFORCE[m, k, j, i], IMFORCE[m, k, j, i] = [float(m) for m in f.readline().split()][1:]

    return REFORCE, IMFORCE, VEL, HEAD, FREQ, XMTN, ZMTN


def read_veres_input(path):
    # TODO: Add documentation

    VMAS, ADDMAS, DAMP, REST, VEL_re7, HEAD_re7, FREQ_re7, XMTN_re7, ZMTN_re7 = read_re7_file(path + '//input.re7')

    REFORCE, IMFORCE, VEL_re8, HEAD_re8, FREQ_re8, XMTN_re8, ZMTN_re8 = read_re8_file(path + '//input.re8')

    A_h = ADDMAS
    B_h = DAMP
    C_h = REST
    F_ex_real = REFORCE
    F_ex_im = IMFORCE

    if VEL_re7.all() == VEL_re8.all() and HEAD_re7.all() == HEAD_re8.all() and FREQ_re7.all() == FREQ_re8.all() and \
            XMTN_re7.all() == XMTN_re8.all() and ZMTN_re7.all() == ZMTN_re8.all():
        VEL = VEL_re7
        HEAD = HEAD_re7
        FREQ = FREQ_re7
        XMTN = XMTN_re7
        ZMTN = ZMTN_re7
    else:
        raise ValueError

    return A_h, B_h, C_h, F_ex_real, F_ex_im, VEL, HEAD, FREQ, XMTN, ZMTN


if __name__ == "__main__":
    path_veres = 'Input files//Veres input files'

    A_h, B_h, C_h, F_ex_real, F_ex_im, VEL, HEAD, FREQ, XMTN, ZMTN = read_veres_input(path_veres)

    print(A_h)

    '''
    REFORCE, IMFORCE, VEL, HEAD, FREQ = read_re8_file('Input files/test_input.re8')

    VMAS, ADDMAS, DAMP, REST, VEL7, HEAD7, FREQ7 = read_re7_file("Input files/test_input.re7")

    print(REFORCE[0, 0, 0, 0])
    '''
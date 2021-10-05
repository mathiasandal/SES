import numpy as np

"""
Reads a *.re7 file created by VERES

A description of a *.re7 file is found in page 114 of the ShipX Vessel Responses (VERES) User's manual.
"""

def read_re7_file(filename):
    f = open(filename, "r")

    for i in range(6):  # skips first six lines
        f.readline()

    #print(isinstance(f.readline(), str))

    #print(f.readline().split())

    #RHOSW, GRAV = f.readline.split()
    #LPP, BREADTH, DRAUGHT = f.readline().split()
    #LCG, VCG = f.readline().split()

    #print([float(i) for i in f.readline().split()])

    # Read in parameters
    RHOSW, GRAV = [float(i) for i in f.readline().split()]
    LPP, BREADTH, DRAUGHT = [float(i) for i in f.readline().split()]
    LCG, VCG = [float(i) for i in f.readline().split()]
    '''
    print(isinstance(RHOSW, float))
    print(isinstance(GRAV, float))
    print()
    print(RHOSW)
    print(GRAV)
    print(LPP)
    print(BREADTH)
    print(DRAUGHT)
    print(LCG)
    print(VCG)
    '''

    # Skip line looking like "-1 FILEVER (=2)" in Veres_Manual.pdf
    f.readline()

    # Read in number of velocities, headings, frequencies and degrees of freedom
    NOVEL, NOHEAD, NOFREQ, NDOF = [int(i) for i in f.readline().split()]

    # Initialize vectors to contain velocities, headings and frequencies
    VEL = np.zeros(NOVEL)
    HEAD = np.zeros(NOHEAD)
    FREQ = np.zeros(NOFREQ)

    # Initialize values to be stores for each velocity
    SINK = np.zeros(NOVEL) # [m] Sink
    TRIM = np.zeros(NOVEL) # [deg] Trim
    XMTN = np.zeros(NOVEL) # [m] X-pos. of the motion coordinate system (relative to Lpp/2)
    ZMTN = np.zeros(NOVEL) # [m] Z-pos. of the mo􀆟on coordinate system (relative to BL)

    # Read in mass matrix
    VMAS = np.zeros([NDOF, NDOF])  # Initialize mass matrix

    for j in range(NDOF):  # Read and converts each line from string to floats
        VMAS[j, :] = [float(i) for i in f.readline().split()]

    # Have one (NDOF x NDOF) matrix for each combinations of velocity, heading and wave frequency
    ADDMAS =    np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Hydrodynamic added mass
    ADDADDMAS = np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Additional added mass
    DAMP =      np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Hydrodynamic damping
    ADDDAMP =   np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Additional damping
    REST =      np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Restoring
    ADDREST =   np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Additional restoring

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
    return VMAS, ADDMAS, DAMP, REST


if __name__ == "__main__":
    VMAS, ADDMAS, DAMP, REST = read_re7_file("Input files/test_input.re7")

    print(np.sqrt(VMAS[4, 4]/VMAS[1, 1]))

from chemex.bases.two_states.iph import R_IXY, DR_IXY, R_IZ, CS, DW, KAB, \
    KBA, W1X


def compute_liouvillian(pb=0.0, kex=0.0, dw=0.0, r_nz=1.5, r_nxy=5.0,
                        dr_nxy=0.0, cs_offset=0.0, w1=0.0):
    """
    Compute the exchange matrix (Liouvillian)

    The function assumes a 2-site (A <-> B) exchanging system.
    The matrix is written in 6x6 cartesian basis, that is {Nx, Ny, Nz}{a,b}.
    Here the thermal equilibrium is assumed to be 0. This is justified
    because of
    the +/- phase cycling of the first 90 degree pulse at the beginning of the
    cpmg block.

    Parameters
    ----------
    pb : float
        Fractional population of state B.
        0.0 for 0%, 1.0 for 100%.
    kex : float
        Exchange rate between state A and B in /s.
    dw : float
        Chemical shift difference between states A and B in rad/s.
    r_nz : float
        Longitudinal relaxation rate of state {a,b} in /s.
    r_nxy : float
        Transverse relaxation rate of state a in /s.
    dr_nxy : float
        Transverse relaxation rate difference between states a and b in /s.
    cs_offset : float
        Offset from the carrier in rad/s.

    Returns
    -------
    out: numpy.matrix
        Liouvillian describing free precession of one
        isolated spin in presence of two-site exchange.

    """

    kab = kex * pb
    kba = kex - kab

    liouvillian = (
        R_IXY * r_nxy +
        DR_IXY * dr_nxy +
        R_IZ * r_nz +
        CS * cs_offset +
        DW * dw +
        w1 * W1X +
        KAB * kab +
        KBA * kba
    )

    return liouvillian
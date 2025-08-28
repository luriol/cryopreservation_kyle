def pt1000_lookup(R_measured):
    """
    Given resistance R in ohms, return the corresponding temperature in °C
    for a Pt1000 using the inverse of the Callendar–Van Dusen equation.
    Returns -999.0 on any error or out-of-range input.
    """
    # Constants for Pt1000 (ITS-90)
    R0 = 1000.0
    A = 3.9083e-3
    B = -5.775e-7
    C_neg = -4.183e-12  # Only used for T < 0 °C

    # Quick input guards
    if R_measured is None or not np.isfinite(R_measured):
        return -999.0

    def R_of_T(T):
        C = C_neg if T < 0 else 0.0
        return R0 * (1 + A*T + B*T*T + C*(T - 100.0)*T**3)

    # Physical validity check for the given bracket
    Tmin, Tmax = -200.0, 850.0
    Rmin, Rmax = R_of_T(Tmin), R_of_T(Tmax)
    lo, hi = (min(Rmin, Rmax), max(Rmin, Rmax))
    if not (lo <= R_measured <= hi):
        return -999.0

    try:
        # Use analytic inverse for T >= 0: R = R0 * (1 + A T + B T^2)
        if R_measured >= R0:
            # Solve B T^2 + A T + (1 - R/R0) = 0
            c = 1.0 - (R_measured / R0)
            disc = A*A - 4.0*B*c
            if disc < 0.0:
                return -999.0
            # Pick the physically meaningful root (positive T, B < 0)
            T1 = (-A + np.sqrt(disc)) / (2.0*B)
            T2 = (-A - np.sqrt(disc)) / (2.0*B)
            # Choose the root ≥ 0 and within [0, Tmax]
            candidates = [t for t in (T1, T2) if np.isfinite(t) and (0.0 <= t <= Tmax)]
            if candidates:
                return float(candidates[0])
            # Fallback to numeric if analytic selection failed
            sol = root_scalar(lambda T: R_of_T(T) - R_measured,
                              bracket=[0.0, Tmax], method='brentq')
            return float(sol.root) if sol.converged else -999.0

        # For T < 0, use full CVD with C term and brentq in [-200, 0]
        sol = root_scalar(lambda T: R_of_T(T) - R_measured,
                          bracket=[Tmin, 0.0], method='brentq')
        return float(sol.root) if sol.converged else -999.0

    except Exception:
        return -999.0

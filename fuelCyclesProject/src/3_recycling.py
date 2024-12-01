from onceThrough_1 import swu

# %% Compute the natural uranium requirement with recycling for an assumed burn-up, enrichment, etc.

# TODO: Essentially the same as before.
# # Coefficient matrix A
# A = np.array([[1, -1, 0],
#               [w_n, -w_t, w_r],
#               [0, 0, 18.3]])
#
# # Right-hand side vector b
# b = np.array([1, w_e, 0])
#
# # Solve the system A x = b
# x = solve(A, b)


# %% Create or use existing SWU function.
htgr_recycling_swu = swu(feed_mass=1, fuel_mass=1, tails_mass=1, feed_assay=1, fuel_assay=1, tails_assay=1)

# %% Input the fuel cost from the website as a constant with recycling.

# %% Compute the masses of the waste stream with recycling.

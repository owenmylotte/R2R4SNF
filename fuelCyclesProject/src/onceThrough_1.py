import numpy as np
from Cython import struct
from scipy.linalg import solve


# %% Compute the natural uranium requirement for an assumed burn-up, enrichment, etc.

class EnrichmentParameters(struct):
    """
    Holds a specific enrichment configuration for a given reactor core.
    """

    def __init__(self, fuel_assay, tails_assay, feed_assay):
        self.fuel_assay = fuel_assay  # Enrichment
        self.tails_assay = tails_assay  # Tails assay
        self.feed_assay = feed_assay  # Natural uranium concentration


# TODO: Need to assume values for the enrichment based on the OpenMC inputs for the reactor materials.
htgr_enrichment_parameters = EnrichmentParameters(fuel_assay=20.0, tails_assay=0.20, feed_assay=0.711)


def input_mass_requirements(enrichment_parameters, mass_fuel):
    """
    Solve for the feed, tails, and fuel requirements given a set mass of required fuel and enrichment parameters.
    """
    # Form the matrix of relationship coefficients.
    A = np.array([[1, -1],
                  [enrichment_parameters.feed_assay, -enrichment_parameters.tails_assay]])

    # Form the right hand side vector.
    b = np.array([1, enrichment_parameters.fuel_assay])

    # Solve the system A x = b for the masses of the feed, fuel, and tails.
    x = solve(A, b)

    return 1 # TODO: This function is incomplete.


# %% Create or use existing SWU function.

def swu(feed_mass, fuel_mass, tails_mass, feed_assay, fuel_assay, tails_assay):
    """
    Outputs the separative work unit (SWU) required for a given enrichment configuration and set of masses.
    """

    # TODO: The input parameters for feed and tails are not required if the input_mass_requirements function is called inside here.
    def swu_value_function(x):
        return (2 * x - 1) * np.log(x / (1 - x))

    return fuel_mass * swu_value_function(fuel_assay) + tails_mass * swu_value_function(
        tails_assay) - feed_mass * swu_value_function(feed_assay)


# Use the function.
htgr_swu = swu(feed_mass=1, fuel_mass=1, tails_mass=1,
               feed_assay=htgr_enrichment_parameters.feed_assay,
               fuel_assay=htgr_enrichment_parameters.fuel_assay,
               tails_assay=htgr_enrichment_parameters.tails_assay)

# %% Input the fuel cost from the website as a constant.

fuel_cost = 0.0  # TODO: Input the fuel cost from the website.

# %% Compute the masses of the waste stream.
# TODO: depleted uranium, spent fuel, effluents and others

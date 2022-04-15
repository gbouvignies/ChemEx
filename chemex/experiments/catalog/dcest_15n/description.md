# 15N Pure In-phase D-CEST

Analyzes chemical exchange in the presence of 1H composite decoupling during
the D-CEST block. This keeps the spin system purely in-phase throughout, and
is calculated using the (3n)×(3n), single-spin matrix, where n is the number
of states:

    { Ix(a), Iy(a), Iz(a),
      Ix(b), Iy(b), Iz(b), ... }

## References

  - Yuwen, Kay and Bouvignies. ChemPhysChem (2018) 19:1707-1710
  - Yuwen, Bah, Brady, Ferrage, Bouvignies and Kay. J Phys Chem B (2018) 122:11206-11217

## Note

A sample configuration file for this module is available using the command:

    $ chemex config dcest_15n

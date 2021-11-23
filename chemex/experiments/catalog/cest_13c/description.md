# 13C Pure In-phase CEST

Analyzes chemical exchange in the presence of 1H composite decoupling during
the CEST block. This keeps the spin system purely in-phase throughout, and is
calculated using the (3n)×(3n), single-spin matrix, where n is the number of
states:

    { Ix(a), Iy(a), Iz(a),
      Ix(b), Iy(b), Iz(b), ... }

## References

  - Vallurupalli, Bouvignies, and Kay. ChemBioChem (2014) 14:1709-1713
  - Bouvignies, Vallurupalli and Kay. J Mol Biol (2014) 426:763-774
  - Vallurupalli and Kay. Angew Chem Int Ed (2013) 52:4156-4159
  - Hansen, Bouvignies and Kay. J Biomol NMR (2013) 55:279-289
  - Bouvignies and Kay. J Biomol NMR (2012) 53:303-310
  - Rennella, Huang, Velyvis and Kay. J Biomol NMR (2015) 63:187-199


## Note

A sample configuration file for this module is available using the command:

    $ chemex config cest_13c

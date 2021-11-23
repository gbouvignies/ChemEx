# 15N Pure In-phase CEST

Analyzes chemical exchange in the presence of 1H composite decoupling during
the CEST block. This keeps the spin system purely in-phase throughout, and is
calculated using the (3n)×(3n), single-spin matrix, where n is the number of
states:

    { Ix(a), Iy(a), Iz(a),
      Ix(b), Iy(b), Iz(b), ... }

## Reference

  - Vallurupalli, Bouvignies and Kay. *J Am Chem Soc* (2012) 134:8148-8161


## Note

A sample configuration file for this module is available using the command:

    $ chemex config cest_15n

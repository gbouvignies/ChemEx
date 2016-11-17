"""Module exchange_model contains code for setting up different three-site
exchange models.
"""


def update_params(params=None,
                  map_names=None,
                  model=None,
                  temperature=None,
                  p_total=None,
                  l_total=None):
    """Update the experimental and fitting parameters depending on the model."""

    if model not in {'pb_kex'}:
        print("Warning: The \'model\' option should be either \'3st.pb_kex\'.")
        print("\nSet it to the default model: \'3st.pb_kex\'.")
        model = 'pb_kex'

"""TODO: module docstring."""

def update_params(params=None, map_names=None, model=None, temperature=None, p_total=None, l_total=None):
    """TODO: function docstring."""
    if model not in {'pb_kex'}:
        print('Warning: The \'model\' option should be either \'3st.pb_kex\'. Set to \'pb_kex\'')
        model = 'pb_kex'

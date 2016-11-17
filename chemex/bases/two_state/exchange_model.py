"""Module exchange_model contains code for setting up different two-site
exchange models.
"""

from scipy import constants as cst

from chemex import parameters


def update_params(params=None,
                  map_names=None,
                  model=None,
                  temperature=None,
                  p_total=None,
                  l_total=None):
    """Update the experimental and fitting parameters depending on the model.

    :param params:
    :param map_names:
    :param model:
    :param temperature:
    :param p_total:
    :param l_total:
    :return:
    """
    if model not in {'2st.pb_kex', '2st.eyring', '2st.kab_kd', '2st.kon_koff'}:
        print("Warning: The \'model\' option should be either \'2st.pb_kex\'")
        print(", \'2st.eyring\', \'2st.kab_kd\', or \'2st.kon_koff\'.")
        print("\nSet it to the default model: \'2st.pb_kex\'.")
        model = '2st.pb_kex'

    if model == '2st.pb_kex':
        pass

    elif model == '2st.eyring':
        map_names.update({
            'dh_b': parameters.ParameterName('dh_b').to_full_name(),
            'ds_b': parameters.ParameterName('ds_b').to_full_name(),
            'dh_ab': parameters.ParameterName('dh_ab').to_full_name(),
            'ds_ab': parameters.ParameterName('ds_ab').to_full_name(),
        })

        t_kelvin = temperature + 273.15
        kbt_h = cst.k * t_kelvin / cst.h
        rt = cst.R * t_kelvin

        pb = ('(1.0 / (1.0 + exp(({dh_b} - {t_kelvin} * {ds_b}) / {rt})))'.format(
            kbt_h=kbt_h, t_kelvin=t_kelvin, rt=rt, **map_names))

        kex_ab = ('{kbt_h} * exp(-({dh_ab} - {t_kelvin} * {ds_ab}) / {rt}) * '
                  '(1.0 + exp(({dh_b} - {t_kelvin} * {ds_b}) / {rt}))'.format(
                      kbt_h=kbt_h, t_kelvin=t_kelvin, rt=rt, **map_names))

        params.add_many(  # Name, Value, Vary, Min, Max, Expr
            (map_names['dh_b'], 6.5e+03, True, None, None, None),
            (map_names['ds_b'], 0.0, False, None, None, None),
            (map_names['dh_ab'], 6.5e+04, True, None, None, None),
            (map_names['ds_ab'], 0.0, False, None, None, None),
            (map_names['pb'], 0.0, None, 0.0, 1.0, pb),
            (map_names['kex_ab'], 0.0, None, 0.0, None, kex_ab), )

    elif model == '2st.kab_kd':
        kwargs1 = {'temperature': temperature}
        kwargs2 = {'temperature': temperature, 'p_total': p_total, 'l_total': l_total}

        map_names.update({
            'kab': parameters.ParameterName('kab', **kwargs1).to_full_name(),
            'kd_ab': parameters.ParameterName('kd_ab', **kwargs1).to_full_name(),
            'l_free': parameters.ParameterName('l_free', **kwargs2).to_full_name(),
        })

        l_free = (
            '0.5 * ({l_total} - {p_total} - {kd_ab} + '
            '(({l_total} - {p_total} - {kd_ab}) ** 2 + 4.0 * {kd_ab} * {l_total})**0.5)'.format(
                p_total=p_total, l_total=l_total, **map_names))

        kex = ('({l_free} + {kd_ab}) * {kab}'.format(p_total=p_total, l_total=l_total, **map_names))

        pb = ('({l_total} - {l_free}) / {p_total}'.format(
            l_total=l_total, p_total=p_total, **map_names))

        params.add_many(  # Name, Value, Vary, Min, Max, Expr
            (map_names['kd_ab'], 100.0, True, 0.0, None, None),
            (map_names['kab'], 1.0e9, True, 0.0, None, None),
            (map_names['l_free'], 0.0, None, 0.0, None, l_free),
            (map_names['pb'], 5.0e-02, None, 0.0, 1.0, pb),
            (map_names['kex_ab'], 2.0e+02, None, 0.0, None, kex), )

    elif model == '2st.kon_koff':
        kwargs1 = {'temperature': temperature}
        kwargs2 = {'temperature': temperature, 'p_total': p_total, 'l_total': l_total}

        map_names.update({
            'kon': parameters.ParameterName('kon', **kwargs1).to_full_name(),
            'koff': parameters.ParameterName('koff', **kwargs1).to_full_name(),
            'kd': parameters.ParameterName('kd', **kwargs1).to_full_name(),
            'l_free': parameters.ParameterName('l_free', **kwargs2).to_full_name(),
        })

        kd = ('{koff} / {kon}'.format(**map_names))

        l_free = ('0.5 * ({l_total} - {p_total} - {kd} + '
                  '(({l_total} - {p_total} - {kd}) ** 2 + 4.0 * {kd} * {l_total})**0.5)'.format(
                      p_total=p_total, l_total=l_total, **map_names))

        kex = ('{l_free} * {kon} + {koff}'.format(p_total=p_total, l_total=l_total, **map_names))

        pb = ('({l_total} - {l_free}) / {p_total}'.format(
            l_total=l_total, p_total=p_total, **map_names))

        params.add_many(  # Name, Value, Vary, Min, Max, Expr
            (map_names['kon'], 1.0e9, True, 0.0, None, None),
            (map_names['koff'], 1.0e2, True, 0.0, None, None),
            (map_names['kd'], 0.0, None, 0.0, None, kd),
            (map_names['l_free'], 0.0, None, 0.0, None, l_free),
            (map_names['pb'], 5.0e-02, None, 0.0, 1.0, pb),
            (map_names['kex_ab'], 2.0e+02, None, 0.0, None, kex), )

    return map_names, params

"""Module exchange_model contains code for setting up different two-site
exchange models."""

import lmfit
from scipy import constants as cst

from chemex import parameters


def create_exchange_params(model=None, temperature=None, p_total=None, l_total=None):
    """Update the experimental and fitting parameters depending on the model.

    :param params:
    :param map_names:
    :param model:
    :param temperature:
    :param p_total:
    :param l_total:
    :return:

    """

    map_names = {}
    params = lmfit.Parameters()

    if model not in {'2st.pb_kex', '2st.eyring', '2st.binding'}:
        print("Warning: The \'model\' option should be either \'2st.pb_kex\'")
        print(", \'2st.eyring\', \'2st.kab_kd\', or \'2st.kon_koff\'.")
        print("\nSet it to the default model: \'2st.pb_kex\'.")
        model = '2st.pb_kex'

    kwargs1 = {'temperature': temperature, 'p_total': p_total, 'l_total': l_total}

    map_names.update({
        'kab': parameters.ParameterName('kab', **kwargs1).to_full_name(),
        'kba': parameters.ParameterName('kba', **kwargs1).to_full_name(),
        'pa': parameters.ParameterName('pa', **kwargs1).to_full_name(),
        'pb': parameters.ParameterName('pb', **kwargs1).to_full_name(),
        'kex_ab': parameters.ParameterName('kex_ab', **kwargs1).to_full_name(),
    })

    if model == '2st.pb_kex':

        params.add_many(
            (map_names['pb'], 0.05, True, 0.0, 1.0, None),
            (map_names['kex_ab'], 200.0, True, 0.0, None, None), )

        pa = '1.0 - {pb}'.format(**map_names)

        params.add_many((map_names['pa'], 0.95, True, 0.0, 1.0, pa), )

        kab = '{kex_ab} * {pb}'.format(**map_names)
        kba = '{kex_ab} * {pa}'.format(**map_names)

        params.add_many(
            (map_names['kab'], 10.0, True, 0.0, None, kab),
            (map_names['kba'], 190.0, True, 0.0, None, kba), )

    elif model == '2st.eyring':
        map_names.update({
            'dh_b': parameters.ParameterName('dh_b').to_full_name(),
            'ds_b': parameters.ParameterName('ds_b').to_full_name(),
            'dh_ab': parameters.ParameterName('dh_ab').to_full_name(),
            'ds_ab': parameters.ParameterName('ds_ab').to_full_name(),
        })

        params.add_many(
            (map_names['dh_b'], 6.5e+03, True, None, None, None),
            (map_names['ds_b'], 0.0, False, None, None, None),
            (map_names['dh_ab'], 6.5e+04, True, None, None, None),
            (map_names['ds_ab'], 0.0, False, None, None, None), )

        t_kelvin = temperature + 273.15
        kbt_h = cst.k * t_kelvin / cst.h
        rt = cst.R * t_kelvin

        kab = ('{kbt_h} * exp(-({dh_ab} - {t_kelvin} * {ds_ab}) / {rt})'.format(
            kbt_h=kbt_h, t_kelvin=t_kelvin, rt=rt, **map_names))

        kba = (
            '{kbt_h} * exp(-(({dh_ab} - {dh_b}) - {t_kelvin} * ({ds_ab} - {ds_b})) / {rt})'.format(
                kbt_h=kbt_h, t_kelvin=t_kelvin, rt=rt, **map_names))

        params.add_many(
            (map_names['kab'], 10.0, True, 0.0, None, kab),
            (map_names['kba'], 190.0, True, 0.0, None, kba), )

        pa = '{kba} / ({kab} + {kba})'.format(**map_names)
        pb = '{kab} / ({kab} + {kba})'.format(**map_names)
        kex_ab = '{kab} + {kba}'.format(**map_names)

        params.add_many(
            (map_names['pa'], 0.0, None, 0.0, 1.0, pa),
            (map_names['pb'], 0.0, None, 0.0, 1.0, pb),
            (map_names['kex_ab'], 0.0, None, 0.0, None, kex_ab), )

    elif model == '2st.binding':
        kwargs1 = {'temperature': temperature}
        kwargs2 = {'temperature': temperature, 'p_total': p_total, 'l_total': l_total}

        map_names.update({
            'kab': parameters.ParameterName('kab', **kwargs1).to_full_name(),
            'kba': parameters.ParameterName('kba', **kwargs1).to_full_name(),
            'kd_ab': parameters.ParameterName('kd_ab', **kwargs1).to_full_name(),
            'l_free': parameters.ParameterName('l_free', **kwargs2).to_full_name(),
        })

        params.add_many(
            (map_names['kab'], 1.0e9, True, 0.0, None, None),
            (map_names['kd_ab'], 100.0, True, 0.0, None, None), )

        kba = '{kd_ab} * {kab}'.format(**map_names)

        l_free = ('0.5 * ({delta} - {kd_ab} + '
                  '(({delta} - {kd_ab}) ** 2 + 4.0 * {kd_ab} * {l_total})**0.5)'.format(
                      delta=l_total - p_total, l_total=l_total, **map_names))

        params.add_many(
            (map_names['kba'], 10.0, None, 0.0, None, kba),
            (map_names['l_free'], 0.0, None, 0.0, None, l_free), )

        pa = '{kba} / ({kab} + {kba})'.format(**map_names)
        pb = '{kab} / ({kab} + {kba})'.format(**map_names)
        kex_ab = '{kab} + {kba}'.format(**map_names)

        params.add_many(
            (map_names['pa'], 0.0, None, 0.0, 1.0, pa),
            (map_names['pb'], 0.0, None, 0.0, 1.0, pb),
            (map_names['kex_ab'], 0.0, None, 0.0, None, kex_ab), )

    return map_names, params

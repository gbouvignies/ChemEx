"""Module exchange_model contains code for setting up different three-site
exchange models."""

import lmfit
from scipy import constants as cst

from chemex import parameters


def create_exchange_params(model=None, temperature=None, p_total=None, l_total=None):
    """Update the experimental and fitting parameters depending on the
    model."""

    map_names = {}
    params = lmfit.Parameters()

    if model not in {'3st.pb_kex', '3st.eyring'}:

        print("Warning: The \'model\' option should be either \'3st.pb_kex\'")
        print(", \'3st.eyring\'.")
        print("\nSet it to the default model: \'3st.pb_kex\'.")

        model = '3st.pb_kex'

    kws = {'temperature': temperature, 'p_total': p_total, 'l_total': l_total}

    map_names.update({
        'kab': parameters.ParameterName('kab', **kws).to_full_name(),
        'kba': parameters.ParameterName('kba', **kws).to_full_name(),
        'kac': parameters.ParameterName('kac', **kws).to_full_name(),
        'kca': parameters.ParameterName('kca', **kws).to_full_name(),
        'kbc': parameters.ParameterName('kbc', **kws).to_full_name(),
        'kcb': parameters.ParameterName('kcb', **kws).to_full_name(),
        'pa': parameters.ParameterName('pa', **kws).to_full_name(),
        'pb': parameters.ParameterName('pb', **kws).to_full_name(),
        'pc': parameters.ParameterName('pc', **kws).to_full_name(),
        'kex_ab': parameters.ParameterName('kex_ab', **kws).to_full_name(),
        'kex_ac': parameters.ParameterName('kex_ac', **kws).to_full_name(),
        'kex_bc': parameters.ParameterName('kex_bc', **kws).to_full_name(),
    })

    if model == '3st.pb_kex':

        params.add_many(
            (map_names['pb'], 0.05, True, 0.0, 1.0, None),
            (map_names['pc'], 0.05, True, 0.0, 1.0, None),
            (map_names['kex_ab'], 200.0, True, 0.0, None, None),
            (map_names['kex_ac'], 0.0, False, 0.0, None, None),
            (map_names['kex_bc'], 200.0, True, 0.0, None, None), )

        pa = '1.0 - {pb} - {pc}'.format(**map_names)

        params.add_many((map_names['pa'], 0.95, True, 0.0, 1.0, pa), )

        kab = '{kex_ab} * {pb} / ({pa} + {pb}) if {pa} + {pb} else 0.0'.format(**map_names)
        kba = '{kex_ab} * {pa} / ({pa} + {pb}) if {pa} + {pb} else 0.0'.format(**map_names)
        kac = '{kex_ac} * {pc} / ({pa} + {pc}) if {pa} + {pc} else 0.0'.format(**map_names)
        kca = '{kex_ac} * {pa} / ({pa} + {pc}) if {pa} + {pc} else 0.0'.format(**map_names)
        kbc = '{kex_bc} * {pc} / ({pb} + {pc}) if {pb} + {pc} else 0.0'.format(**map_names)
        kcb = '{kex_bc} * {pb} / ({pb} + {pc}) if {pb} + {pc} else 0.0'.format(**map_names)

        params.add_many(
            (map_names['kab'], 10.0, True, 0.0, None, kab),
            (map_names['kba'], 190.0, True, 0.0, None, kba),
            (map_names['kac'], 0.0, True, 0.0, None, kac),
            (map_names['kca'], 0.0, True, 0.0, None, kca),
            (map_names['kbc'], 0.0, True, 0.0, None, kbc),
            (map_names['kcb'], 0.0, True, 0.0, None, kcb), )

    elif model == '3st.eyring':

        map_names.update({
            'dh_b': parameters.ParameterName('dh_b').to_full_name(),
            'ds_b': parameters.ParameterName('ds_b').to_full_name(),
            'dh_c': parameters.ParameterName('dh_c').to_full_name(),
            'ds_c': parameters.ParameterName('ds_c').to_full_name(),
            'dh_ab': parameters.ParameterName('dh_ab').to_full_name(),
            'ds_ab': parameters.ParameterName('ds_ab').to_full_name(),
            'dh_ac': parameters.ParameterName('dh_ac').to_full_name(),
            'ds_ac': parameters.ParameterName('ds_ac').to_full_name(),
            'dh_bc': parameters.ParameterName('dh_bc').to_full_name(),
            'ds_bc': parameters.ParameterName('ds_bc').to_full_name(),
        })

        params.add_many(
            (map_names['dh_b'], 6.5e+03, True, None, None, None),
            (map_names['ds_b'], 0.0, False, None, None, None),
            (map_names['dh_c'], 6.5e+03, True, None, None, None),
            (map_names['ds_c'], 0.0, False, None, None, None),
            (map_names['dh_ab'], 6.5e+04, True, None, None, None),
            (map_names['ds_ab'], 0.0, False, None, None, None),
            (map_names['dh_ac'], 6.5e+09, False, None, None, None),
            (map_names['ds_ac'], 0.0, False, None, None, None),
            (map_names['dh_bc'], 6.5e+04, True, None, None, None),
            (map_names['ds_bc'], 0.0, False, None, None, None), )

        t_kelvin = temperature + 273.15
        kbt_h = cst.k * t_kelvin / cst.h
        rt = cst.R * t_kelvin

        kab = ('{kbt_h} * exp(-({dh_ab} - {t_kelvin} * {ds_ab}) / {rt})'.format(
            kbt_h=kbt_h, t_kelvin=t_kelvin, rt=rt, **map_names))

        kba = (
            '{kbt_h} * exp(-(({dh_ab} - {dh_b}) - {t_kelvin} * ({ds_ab} - {ds_b})) / {rt})'.format(
                kbt_h=kbt_h, t_kelvin=t_kelvin, rt=rt, **map_names))

        kac = ('{kbt_h} * exp(-({dh_ac} - {t_kelvin} * {ds_ac}) / {rt})'.format(
            kbt_h=kbt_h, t_kelvin=t_kelvin, rt=rt, **map_names))

        kca = (
            '{kbt_h} * exp(-(({dh_ac} - {dh_c}) - {t_kelvin} * ({ds_ac} - {ds_c})) / {rt})'.format(
                kbt_h=kbt_h, t_kelvin=t_kelvin, rt=rt, **map_names))

        kbc = (
            '{kbt_h} * exp(-(({dh_bc} - {dh_b}) - {t_kelvin} * ({ds_bc} - {ds_b})) / {rt})'.format(
                kbt_h=kbt_h, t_kelvin=t_kelvin, rt=rt, **map_names))

        kcb = (
            '{kbt_h} * exp(-(({dh_bc} - {dh_c}) - {t_kelvin} * ({ds_bc} - {ds_c})) / {rt})'.format(
                kbt_h=kbt_h, t_kelvin=t_kelvin, rt=rt, **map_names))

        params.add_many(
            (map_names['kab'], 10.0, True, 0.0, None, kab),
            (map_names['kba'], 190.0, True, 0.0, None, kba),
            (map_names['kac'], 0.0, True, 0.0, None, kac),
            (map_names['kca'], 0.0, True, 0.0, None, kca),
            (map_names['kbc'], 0.0, True, 0.0, None, kbc),
            (map_names['kcb'], 0.0, True, 0.0, None, kcb), )

        pa = '{kba} * {kca} / ({kba} * {kca} + {kab} * {kca} + {kac} * {kba})'.format(**map_names)
        pb = '{kab} * {kcb} / ({kab} * {kcb} + {kba} * {kcb} + {kbc} * {kab})'.format(**map_names)
        pc = '{kbc} * {kac} / ({kbc} * {kac} + {kcb} * {kab} + {kca} * {kbc})'.format(**map_names)
        kex_ab = '{kab} + {kba}'.format(**map_names)
        kex_ac = '{kac} + {kca}'.format(**map_names)
        kex_bc = '{kbc} + {kcb}'.format(**map_names)

        params.add_many(
            (map_names['pa'], 0.0, None, 0.0, 1.0, pa),
            (map_names['pb'], 0.0, None, 0.0, 1.0, pb),
            (map_names['pc'], 0.0, None, 0.0, 1.0, pc),
            (map_names['kex_ab'], 0.0, None, 0.0, None, kex_ab),
            (map_names['kex_ac'], 0.0, None, 0.0, None, kex_ac),
            (map_names['kex_bc'], 0.0, None, 0.0, None, kex_bc), )

    return map_names, params

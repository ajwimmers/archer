from __future__ import print_function
import numpy as np
import archer.utilities.InterpToolbox as intbx
#from importlib import reload
#reload(intbx)


def confidence_to_alpha(confidence_score, archer_channel_type, fxHr, vmax):

    alpha_floor = 0.5

    if archer_channel_type == 'RSCAT':
        m_fit = 1.8
        b_fit = -19.40

    elif archer_channel_type == 'ASCAT':
        m_fit = 2.67
        b_fit = -19.40

    elif archer_channel_type == '89GHz':
        m_fit_lo = 3.28
        b_fit_lo = 2.63
        m_fit_hi = 1.61
        b_fit_hi = 9.58

    elif archer_channel_type == '37GHz':
        m_fit_lo = 3.28
        b_fit_lo = 2.63
        m_fit_hi = 1.61
        b_fit_hi = 9.58

    elif archer_channel_type == '183GHz':
        m_fit_lo = 3.28
        b_fit_lo = 2.63
        m_fit_hi = 1.61
        b_fit_hi = 9.58

    elif archer_channel_type == 'IR':
        m_fit_lo = 9.89
        b_fit_lo = -2.07
        m_fit_hi = 9.26
        b_fit_hi = 1.95

    elif archer_channel_type == 'SWIR':
        m_fit_lo = 8.68
        b_fit_lo = -0.37
        m_fit_hi = 14.2
        b_fit_hi = -0.24

    elif archer_channel_type in ['Vis', 'DNB']:
        m_fit_lo = 14.44
        b_fit_lo = -0.83
        m_fit_hi = 14.64
        b_fit_hi = 3.45


    # Compute alpha_0 (alpha at dt=0) from the above parameters
    if archer_channel_type in ['RSCAT', 'ASCAT']:
        alpha_0 = np.max(m_fit * confidence_score + b_fit, alpha_floor)

    elif archer_channel_type in ['37GHz', '89GHz', '183GHz', 'IR', 'SWIR', 'Vis', 'DNB']:
        alpha_lo0 = np.max([m_fit_lo * confidence_score + b_fit_lo, alpha_floor])
        alpha_hi0 = np.max([m_fit_hi * confidence_score + b_fit_hi, alpha_floor])

        alpha_0 = np.interp(vmax, [0, 60, 85, 300], [alpha_lo0, alpha_lo0, alpha_hi0, alpha_hi0])

    elif archer_channel_type == 'NHCanfx':
        alpha_at_48_kt = 9.5
        alpha_0 = alpha_at_48_kt + 0.05 * (vmax - 48)


    # Compute alpha at dt = fxHr for either analysis/forecast or imagery
    if archer_channel_type == 'NHCanfx':
        c1 = 3.00/13.4/(15**1.5)
        alpha = alpha_0 / (1 + c1 * alpha_0 * (fxHr**1.50) )

    else:
        c1 = 4.00/13.4/(15**1.5)
        alpha = alpha_0 / (1 + c1 * alpha_0 * (fxHr**1.50) )

    return alpha


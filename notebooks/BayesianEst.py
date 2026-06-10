import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, norm

def get_est(df_iCESM, mdl_TEMP, mdl_d18Osw, d18Ob_mu, d18Ob_sigma, d18Osw_mu=-0.55, d18Osw_sigma=0.1, DpH_mu=-0.2, DpH_sigma=0.05):
    # global mean
    dBWT = {}
    dGMST = {}
    dGMSST = {}

    # benthic sites average
    TEMP = {}
    d18Osw = {}

    for tag in ['1.5x', '3x']:
        dBWT[tag] = df_iCESM.loc[tag, 'dBWT']
        dGMST[tag] = df_iCESM.loc[tag, 'dGMST']
        dGMSST[tag] = df_iCESM.loc[tag, 'dGMSST']
        TEMP[tag] = mdl_TEMP[tag]
        d18Osw[tag] = mdl_d18Osw[tag]

    # extrpolate to get 6x
    dBWT['6x'] = dBWT['3x'] + (dBWT['3x'] - dBWT['1.5x'])
    dGMST['6x'] = dGMST['3x'] + (dGMST['3x'] - dGMST['1.5x'])
    dGMSST['6x'] = dGMSST['3x'] + (dGMSST['3x'] - dGMSST['1.5x'])
    TEMP['6x'] = TEMP['3x'] + (TEMP['3x'] - TEMP['1.5x'])
    d18Osw['6x'] = d18Osw['3x'] + (d18Osw['3x'] - d18Osw['1.5x'])

    # interpolation & extrapolation of the iCESM1.3 simualtions
    n = 500 
    a, b = 1.5, 3
    log_vals1 = np.logspace(np.log10(a), np.log10(b), num=n)

    a, b = 3, 6
    log_vals2 = np.logspace(np.log10(a), np.log10(b), num=n)

    CO2_levs = np.append(log_vals1, log_vals2[1:])

    # climate variables

    dGMST_vals = np.linspace(dGMST['1.5x'], dGMST['6x'], len(CO2_levs))
    dGMSST_vals = np.linspace(dGMSST['1.5x'], dGMSST['6x'], len(CO2_levs))
    dBWT_vals = np.linspace(dBWT['1.5x'], dBWT['6x'], len(CO2_levs))
    TEMP_vals = np.linspace(TEMP['1.5x'], TEMP['6x'], len(CO2_levs))
    d18Osw_vals = np.linspace(d18Osw['1.5x'], d18Osw['6x'], len(CO2_levs))

    def d18Ob_M14(T, d18Osw, DpH):
        d18Ob = (d18Osw-0.27) + (-0.245*T+0.0011*T*T + 3.58) - 1.5435*DpH  # -0.27: VSMOW to VPDB conversion 
        return d18Ob

    d18Ob_vals = d18Ob_M14(TEMP_vals, d18Osw_vals, 0)  # internal d18Osw, no pH correction

    # S = d18Obase - 1.5435*DpH
    k = 1.5435
    S_mu = d18Osw_mu - k*DpH_mu
    S_sigma = np.sqrt(d18Osw_sigma**2 + (k*DpH_sigma)**2)
    centers = S_mu + d18Ob_vals
    area = norm.cdf(d18Ob_mu + d18Ob_sigma, centers, S_sigma) \
         - norm.cdf(d18Ob_mu - d18Ob_sigma, centers, S_sigma)
    area_pdf = area / area.sum()

    d18Ob_x = np.linspace(0, 3, n)
    density_array = norm.pdf(d18Ob_x[None, :], centers[:, None], S_sigma)

    x = np.log(CO2_levs)
    kde = gaussian_kde(x, weights=area_pdf)
    pdf = kde(x)
    df_CO2pdf = pd.DataFrame({'CO2': CO2_levs, 'density': pdf, 'dBWT': dBWT_vals, 'dGMST': dGMST_vals, 'dGMSST': dGMSST_vals})

    res = {
        'density_array': density_array,
        'd18Ob_x': d18Ob_x,
        'df_CO2pdf': df_CO2pdf,
    }
    return res
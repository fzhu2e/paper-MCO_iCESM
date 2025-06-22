import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.integrate as integrate
from scipy.integrate import simpson
from scipy.stats import gaussian_kde
import pymc as pm

def get_est(d18Osw_mu=-0.55, d18Osw_sigma=0.1, DpH_mu=-0.2, DpH_sigma=0.05):
    df_sites_MCO = pd.read_csv('../data/df_d18Ob_sites_MCO.csv')
    df_sites = pd.read_csv('../data/df_d18Ob_sites_model_vals.csv')
    mask_pid = df_sites['pid'].isin(df_sites_MCO['pid'])
    df_sites = df_sites[mask_pid].set_index('pid')

    df_eq_gmst = pd.read_csv('../data/df_eq_gmst.csv', index_col=0)
    df_eq_bwt = pd.read_csv('../data/df_eq_bwt.csv', index_col=0)
    df_eq_gmsst = pd.read_csv('../data/df_eq_gmsst.csv', index_col=0)

    # global mean
    dBWT = {}
    dGMST = {}
    dGMSST = {}

    # benthic sites average
    TEMP = {}
    d18Osw = {}

    for tag in ['PI', '1.5x', '3x']:
        dBWT[tag] = df_eq_bwt.loc[tag, 'dBWT']
        dGMST[tag] = df_eq_gmst.loc[tag, 'dGMST']
        dGMSST[tag] = df_eq_gmsst.loc[tag, 'dGMSST']
        TEMP[tag] = df_sites[f'TEMP ({tag})'].mean()
        d18Osw[tag] = df_sites[f'd18Osw ({tag})'].mean()

    # extrpolate to get 6x
    dBWT['6x'] = dBWT['3x'] + (dBWT['3x'] - dBWT['1.5x'])
    dGMST['6x'] = dGMST['3x'] + (dGMST['3x'] - dGMST['1.5x'])
    dGMSST['6x'] = dGMSST['3x'] + (dGMSST['3x'] - dGMSST['1.5x'])
    TEMP['6x'] = TEMP['3x'] + (TEMP['3x'] - TEMP['1.5x'])
    d18Osw['6x'] = d18Osw['3x'] + (d18Osw['3x'] - d18Osw['1.5x'])


    # observational uncertainty range
    df_proxy_Bd18O = pd.read_csv('../data/df_proxy_d18Ob.csv', index_col=0)
    upper_bd = df_proxy_Bd18O.loc['Dd18Ob', 'upper']
    lower_bd = df_proxy_Bd18O.loc['Dd18Ob', 'lower']


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
        d18Ob = (d18Osw-0.27) + (-0.245*T+0.0011*T*T + 3.58) - 1.5435*DpH
        return d18Ob

    d18Ob_vals = d18Ob_M14(TEMP_vals, d18Osw_vals, 0)  # internal d18Osw, no pH correction

    PI_d18Ob = df_sites['d18Ob (PI)'].mean()
    d18Ob_vals -= PI_d18Ob

    # Lear et al. (2015)
    with pm.Model() as model:
        d18Osw = pm.Normal.dist(mu=d18Osw_mu, sigma=d18Osw_sigma)

    # S = d18Obase - 1.5435*DpH
    k = 1.5435
    S_mu = d18Osw_mu - k*DpH_mu
    S_sigma = np.sqrt(d18Osw_sigma**2 + (k*DpH_sigma)**2)
    area = d18Ob_vals.copy()

    S_pdf_list = []
    for i in range(len(d18Ob_vals)):
        def S_pdf(x):
            return (1 / (S_sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - (S_mu + d18Ob_vals[i])) / S_sigma) ** 2)

        S_pdf_list.append(S_pdf)
        area[i], _ = integrate.quad(S_pdf, lower_bd, upper_bd)

    area_pdf = area / area.sum()

    x = np.linspace(-3, 0, n)

    density_array = np.ndarray((len(CO2_levs), len(x)))
    for i in range(len(CO2_levs)):
        density_array[i] = S_pdf_list[i](x)

    # CO2 --------------------------------------
    x = np.log(CO2_levs)
    kde = gaussian_kde(x, weights=area_pdf)
    pdf = kde(x)
    df_CO2pdf = pd.DataFrame({'CO2': CO2_levs, 'density': pdf, 'dBWT': dBWT_vals, 'dGMST': dGMST_vals, 'dGMSST': dGMSST_vals})


    res = {}
    # dBWT --------------------------------------
    x = np.linspace(0, 15, 151)
    data = df_CO2pdf['dBWT'].values
    kde_iCESM = gaussian_kde(data, weights=df_CO2pdf['density'])
    pdf_iCESM = kde_iCESM(x)
    prior = pdf_iCESM
    mle = x[np.argmax(prior)]

    # Mean
    mean = simpson(y=x * pdf_iCESM, x=x)
    # Variance
    variance = simpson(y=(x - mean)**2 * pdf_iCESM, x=x)
    # Standard Deviation
    std_dev = np.sqrt(variance)

    res['dBWT'] = {
        'mean': mean,
        'std_dev': std_dev,
        'mle': mle,
    }

    # dGMST --------------------------------------
    x = np.linspace(0, 15, 151)
    data = df_CO2pdf['dGMST'].values
    kde_iCESM = gaussian_kde(data, weights=df_CO2pdf['density'])
    pdf_iCESM = kde_iCESM(x)
    prior = pdf_iCESM
    mle = x[np.argmax(prior)]

    # Mean
    mean = simpson(y=x * pdf_iCESM, x=x)
    # Variance
    variance = simpson(y=(x - mean)**2 * pdf_iCESM, x=x)
    # Standard Deviation
    std_dev = np.sqrt(variance)

    res['dGMST'] = {
        'mean': mean,
        'std_dev': std_dev,
        'mle': mle,
        'df_CO2pdf': df_CO2pdf,
    }

    return res
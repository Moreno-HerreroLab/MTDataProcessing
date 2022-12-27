import numpy as np
from scipy.optimize import curve_fit
from mod import process_data as prc

def fit_peak(bins_center, counts, limits, common, sigmas):
    mask, _, _ = limits
    sigma0, sigmamin, sigmamax = sigmas
    x_to_fit, y_to_fit = bins_center[mask], counts[mask]
    try:
        def Gauss(x, a, mu, sigma):
            return a*np.exp(-(x-mu)**2/(2*sigma**2))
        xpeak, ypeak = prc.find_coord_hist_maxima(x_to_fit, y_to_fit)
        #params = gauss_fit(Gauss, x_to_fit, y_to_fit, limits, x_peak, y_peak)
        p0 = ypeak, xpeak, sigma0   # a0, mu0, sigma0
        params = gauss_fit(Gauss, x_to_fit, y_to_fit, limits, p0, sigmas)
        a, _, sigma = params
        if sigma > 1 or a < 0.5: 
            print('Not a peak, too small')
            return None
        else: 
            percent = plot_fitting(Gauss, params, bins_center, common)
            return percent
    except: 
        print('Unable to fit.')
        return None


#def gauss_fit(Gauss, x_to_fit, y_to_fit, limits, x_peak, y_peak):
def gauss_fit(Gauss, x_to_fit, y_to_fit, limits, p0, sigmas):
    a0, mu0, sigma0 = p0
    _, sigmamin, sigmamax = sigmas
    _, zmin, zmax = limits
    #sigma = sum(y_to_fit*(x_to_fit-mu)**2)
    min = [a0-0.1, zmin, sigmamin]
    max = [a0+0.1, zmax, sigmamax]
    #params, cov = curve_fit(Gauss, x_to_fit, y_to_fit, p0=[a, mu, sigma], bounds=(min, max))
    params, cov = curve_fit(Gauss, x_to_fit, y_to_fit, p0=p0, bounds=(min, max))
    return params


def plot_fitting(Gauss, params, bins_center, common):
    ax, bw = common
    a, mu, sigma = params
    # plot the gaussian function in the "entire" x range
    xrange = np.linspace(-0.6, 1, 5000)
    x, y = xrange, Gauss(xrange,*params)
    ax.plot(x, y, 'k-', linewidth=1)
    # annotate params    
    ax.annotate(f'{round(mu, 3)} µm', xy=(mu, 3), size=10, ha='center', va='center', rotation=90)
    print(f'a: {round(a, 3)}, mu: {round(mu, 3)}, sigma: {round(sigma, 3)}, ', end = '')
    # find the area under the curve, actually the area for the hist rectangles
    x, y = bins_center, Gauss(bins_center,*params)
    area = sum(y*bw)
    percent = int(area*100) # the histogram is normalized. the sum of all area is 1. (area_peak*100)/1
    print('area %:', percent)
    ax.annotate(f'{percent} %', xy=(mu, 0.5), xycoords='data', size=10, ha='center', va='center')
    return percent


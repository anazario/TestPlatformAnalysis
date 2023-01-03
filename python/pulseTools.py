from scipy import signal
from scipy import interpolate
from scipy.ndimage import filters
from scipy import optimize

import BSpline as bsp
import matplotlib.pyplot as plt
import numpy as np
import math
import ROOT

def cfd(waveform, time, threshold=0.02, base_line=20, fraction=0.3, hysteresis=0.001, smooth=False, kernel_width=5, kernel_sigma=1):
    samples = np.copy(waveform)
    if samples.max()-samples.min() < threshold:
        return None
    if smooth:
        kernel = signal.gaussian(kernel_width, kernel_sigma)
        samples = filters.convolve(samples, kernel)
    samples -= np.mean(samples[:base_line])
    #print(samples)
    if abs(np.max(samples)) < abs(np.min(samples)):
        maximum = abs(samples.min())
        samples /= -maximum
    else:
        maximum = samples.max()
        samples /= maximum
    #print(samples)
    hysteresis = hysteresis/maximum
    above = False
    lockedForHysteresis = False
    numberOfPeaks = 0
    for i,v in enumerate(samples):
        if not above and not lockedForHysteresis and v>fraction:
            firstCrossingIndex = i
            above = True
            lockedForHysteresis = True
        if above and lockedForHysteresis and v>fraction+hysteresis:
            lockedForHysteresis = False
        if above and not lockedForHysteresis and v<fraction:
            numberOfPeaks += 1
            above = False

    if numberOfPeaks==1:
        #print(f'crossing index: {firstCrossingIndex}')
        t1 = time[firstCrossingIndex-1]
        t2 = time[firstCrossingIndex]
        v1 = samples[firstCrossingIndex-1]
        v2 = samples[firstCrossingIndex]
        #print(f't1: {t1} ns, t2: {t2} ns')
        #print(f'v1: {v1}, v2: {v2}')
        return t1 + (t2 - t1)*(fraction-v1)/(v2 - v1)
    else:
        return None

def cfd_aligned_waveform(pulse, size, cfd_fraction, baseline):
    
    cfd_time = pulse.get_cfd_time(fraction=cfd_fraction, baseline=baseline)#, smooth=True)                                                                            
    transformation_time = np.zeros(size)
    transformed_waveform = np.zeros(size)
    
    if cfd_time is not None and abs(cfd_time) < 5:                                                                                                                  
        transformation_time = np.linspace(cfd_time-1.8, cfd_time+20, size)
        transformation_function = interpolate.interp1d(x=pulse.get_time_array(), y=pulse.get_waveform(), kind='cubic')
        transformed_waveform = transformation_function(transformation_time)
    
    return transformation_time, transformed_waveform
        

def get_interpolated_point(pulse, time, frequency):

    size = len(pulse)

    total_sum = 0
    for i in range(size):
        argument = math.pi*(frequency*time-float(i))
        sinc = -999.

        if abs(argument) < 1e-8:
            sinc = 1.0
        else:
            sinc = math.sin(argument)/argument

        total_sum += pulse[i]*sinc

    if math.isnan(total_sum):
        return 0
    else:
        return total_sum

def poly_func(coeff_array, x_val):
    total = 0
    arr_size = len(coeff_array)-1
    for coeff in coeff_array:
        total += coeff*x_val**arr_size
        arr_size -= 1
    return total

def exp_func(x, a, b, c):
    return a * np.exp(-b * x) + c

def segments_fit(X, Y, maxcount):

    X = np.array(X)
    Y = np.array(Y)
    
    xmin = np.min(X)
    xmax = np.max(X)
    
    n = len(X)
    
    AIC_ = float('inf')
    BIC_ = float('inf')
    r_   = None
    
    for count in range(1, maxcount+1):
        
        seg = np.full(count - 1, (xmax - xmin) / count)

        px_init = np.r_[np.r_[xmin, seg].cumsum(), xmax]
        py_init = np.array([Y[np.abs(X - x) < (xmax - xmin) * 0.1].mean() for x in px_init])

        def func(p):
            seg = p[:count - 1]
            py = p[count - 1:]
            px = np.r_[np.r_[xmin, seg].cumsum(), xmax]
            return px, py

        def err(p): # This is RSS / n
            px, py = func(p)
            Y2 = np.interp(X, px, py)
            return np.mean((Y - Y2)**2)

        r = optimize.minimize(err, x0=np.r_[seg, py_init], method='Nelder-Mead')
    
        # Compute AIC/ BIC. 
        AIC = n * np.log10(err(r.x)) + 4 * count
        BIC = n * np.log10(err(r.x)) + 2 * count * np.log(n)
        
        if (BIC < BIC_) and (AIC < AIC_): # Continue adding complexity.
            r_ = r
            AIC_ = AIC
            BIC_ = BIC
        else: # Stop.
            count = count - 1
            break
        
    return func(r_.x) ## Return the last (n-1)

def inv_piecewise_func(x, line_x_arr, line_y_arr):
    length = len(line_x_arr)
    coeff_collection = []

    for i in range(1,length):
        if x > line_x_arr[i-1] and x < line_x_arr[i]:
            pair_x = [line_x_arr[i-1], line_x_arr[i]]
            pair_xidx = [i-1, i]
            break
    
    for i in range(1,length):
        pair_x = [line_x_arr[i-1], line_x_arr[i]]
        pair_y = [line_y_arr[i-1], line_y_arr[i]]
        coeff_collection.append(np.polyfit(pair_x, pair_y, 1))
    print(coeff_collection)
    return 0

def bspline(samples, knot_vector, order=3, precision=1e-8, isClamped=False):

    fit = bsp.BSpline(samples, knot_vector, order=order, precision=precision, isClamped=isClamped)
    spline = fit.spline_function()
    error = fit.get_error()

    return spline, error

def bspline_err(knot_vector, samples):
    knot_vector = list(knot_vector)
    knot_vector.sort()
    fit = bsp.BSpline(samples, knot_vector)
    error = fit.get_error()

    return error
        
def plot_scatter(x_array, y_array, name, xlabel, ylabel, title='', fit=None):
    fig, ax = plt.subplots()
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    plt.grid(linestyle = '--', linewidth = 0.5)
    plt.scatter(x_array, y_array, color='black', s=4)

    if fit is not None:
        #fit_coeff = np.polyfit(x_array, y_array, fit)
        #fit_coeff, fit_cov = curve_fit(exp_func, x_array, y_array, maxfev=len(x_array))
        #print(fit_coeff)
        #x_fit = np.linspace(np.min(x_array), np.max(x_array), 200)
        #y_fit = poly_func(fit_coeff, x_fit)
        x_fit, y_fit = segments_fit(x_array, y_array, fit)
        print(x_fit, y_fit)
        plt.plot(x_fit, y_fit)
        
    plt.savefig(name+'.pdf')

    print(f'Saved plot: {name}.pdf')

def plot_line(x_array, y_array, name, xlabel, ylabel, title, leg_labels=None):
    colors = ['black','red', 'blue', 'green', 'orange']
    x_array = np.array(x_array)
    y_array = np.array(y_array)

    fig, ax = plt.subplots()
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    plt.grid(linestyle = '--', linewidth = 0.5)

    if y_array.ndim == 2:
        for idx, yarr in enumerate(y_array):
            if leg_labels is not None:
                plt.plot(x_array, yarr, color=colors[idx], label=leg_labels[idx])
            else:
                plt.plot(x_array, yarr, color=colors[idx])
    elif y_array.ndim == 1:
        plt.plot(x_array, y_array, color='black')
    else:
        sys.exit('Incorrect input format!')

    if leg_labels is not None:
        plt.legend()
        
    plt.savefig(name+'.pdf')
    print(f'Saved plot: {name}.pdf')

def plot_hist(data, bins, range_tuple, name, xlabel, ylabel, title):
    fig, ax = plt.subplots()
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    plt.grid(linestyle = '--', linewidth = 0.5)

    plt.hist(data, bins=bins, range=range_tuple)
    plt.savefig(name+'.pdf')
    print(f'Saved plot: {name}.pdf')


#ROOT plotting macros

def plot_2D_hist(hist_th2, name, x_label, y_label, setLogz=False):
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    c1 = ROOT.TCanvas("c1","c1",600,800)
    if setLogz:
        c1.SetLogz()
    c1.SetLeftMargin(0.12);
    c1.SetRightMargin(0.18)
    c1.SetGrid()
    hist_th2.GetXaxis().SetTitle(x_label)
    hist_th2.GetXaxis().CenterTitle()
    hist_th2.GetYaxis().SetTitle(y_label)
    hist_th2.GetYaxis().CenterTitle()
    hist_th2.GetZaxis().SetTitle('Events')
    hist_th2.GetZaxis().CenterTitle()
    hist_th2.GetZaxis().SetTitleOffset(1.5)
    hist_th2.Draw("COLZ")
    c1.SaveAs(f'{name}.pdf')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def main():
    fname = option_dict['fname']
    plot_2x2_heatmap(fname)
    plot_single(fname,2)

def format_data(data,channel):
    x = data[:,1].round(2)
    y = data[:,0].round(2)
    ch = data[:,channel+1]
    #ch = [1./c if c > 0 else 0 for c in ch]
    
    dfch = pd.DataFrame.from_dict(np.array([x,y,ch]).T)
    dfch.columns = ['x','y','z']
    dfch['z'] = pd.to_numeric(dfch['z'])
    pvt = dfch.pivot('y','x','z')

    return pvt

def make_heatmap(pvt,axis,channel):
    axs = axis[f'{channel}']
    sns.heatmap(pvt[f'ch{channel}'] ,cmap=cmap_opt,ax=axs,cbar_kws={'label': 'rate [Hz]'})

    axs.set_title(f'ch{channel}')
    axs.set_xlabel('y [mm]')
    axs.set_ylabel('x [mm]')
    #axs.invert_yaxis()

def plot_single(fname,channel):
    data = np.loadtxt(fname)
    fname_tmp = fname.replace('data/','')
    num_row, num_col = data.shape

    pvt = format_data(data,channel)

    fig, ax = plt.subplots(constrained_layout=True)
    sns.heatmap(pvt,cmap=cmap_opt,ax=ax,cbar_kws={'label': 'rate [Hz]'})
    ax.set_xlabel('y [mm]')
    ax.set_ylabel('x [mm]')
    #ax.invert_yaxis()
    fig.suptitle('Platform Scan Beam Profile '+f'ch{channel}')
    
    plt_name = 'plots/beamProfile_'+fname_tmp.replace('.txt','_')+f'ch{channel}.pdf'
    plt.savefig(plt_name)
    print("Created file "+plt_name)

def plot_2x2_heatmap(fname):
    data = np.loadtxt(fname)
    fname_tmp = fname.replace('data/','')
    num_row, num_col = data.shape

    #prepare axes for plotting heatmaps
    ax = dict.fromkeys([f'{i}' for i in range(1,num_col-1)])
    fig, ((ax['4'], ax['3']), (ax['1'], ax['2'])) = plt.subplots(2, 2, constrained_layout=True)

    pvt = {}
    for ch in range(1,num_col-1):
        pvt[f'ch{ch}'] = (format_data(data,ch))
        make_heatmap(pvt,ax,ch)

    fig.suptitle('Platform Scan Beam Profile')
    plt_name = 'plots/beamProfile_'+fname_tmp.replace('.txt','')+'.pdf'
    plt.savefig(plt_name)
    print("Created file "+plt_name)
    
if __name__ == "__main__":
    import optparse

    usage = "usage: %prog -f <filename>"
    parser = optparse.OptionParser(usage, version="%prog 0.1.0")
    parser.add_option("-f", "--file", type="str",
                      dest="fname",
                      help="Specify data input file",
                      default="scan11.txt")
    (options, args) = parser.parse_args()
    option_dict = vars(options)

    #heatmap style option
    cmap_opt = 'afmhot'
    
    main()

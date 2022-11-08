from functionsmodule import *

import pandas_bokeh
from venn import venn
import matplotlib.pyplot as plt


sample = sys.argv[1]
data = sys.argv[2]
path = './output/'
sites = glob.glob(path + sample +'/' + data)

if data == 'filtered_variants.xlsx':
    pontype = ""
elif data == 'filtered_variants_PoN.xlsx':
    pontype = "PoN"
else:
  print("Wrong input")

df = pd.read_excel(sites[0], engine='openpyxl')
df['TOOLS'] = df['Samples'].str.split('VD:', 1).str[1].str.split('_AF', 1).str[0]
df['AF_mean'] = df['Samples'].str.split('AF:', 1).str[1].str.split('_DP', 1).str[0].astype(float)
df['DP_mean'] = df['Samples'].str.split('DP:', 1).str[1].str.split(';', 1).str[0].astype(int)


### GENERATING VENN
dfvenn = df.copy()
dfvenn['SV'] = dfvenn['TOOLS'].str[0:1]
dfvenn['M2'] = dfvenn['TOOLS'].str[1:2]
dfvenn['LF'] = dfvenn['TOOLS'].str[2:3]
dfvenn['VD'] = dfvenn['TOOLS'].str[3:4]
dfvenn['loc'] = dfvenn['Chrom'] + '_' + dfvenn['Position'].astype(str) + '_' + dfvenn['Ref Base'] + '>' + dfvenn['Alt Base']

dfSV = dfvenn.loc[dfvenn['SV'] == '1']
SVlist = list(dfSV['loc'])
dfM2 = dfvenn.loc[dfvenn['M2'] == '1']
M2list = list(dfM2['loc'])
dfLF = dfvenn.loc[dfvenn['LF'] == '1']
LFlist = list(dfLF['loc'])
dfVD = dfvenn.loc[dfvenn['VD'] == '1']
VDlist = list(dfVD['loc'])

SNVcalls = {
    "SiNVICT": {i for i in SVlist},
    "Mutect2": {i for i in M2list},
    "LoFreq": {i for i in LFlist},
    "VarDict": {i for i in VDlist}
}

fig = venn(SNVcalls
     #, fmt="{percentage:.1f}%"
    ).figure
fig.suptitle(sys.argv[1] + "  " + pontype, fontsize=15)
plt.xlabel('Total # of mutations: '+ str(len(dfvenn)), fontsize=12)
fig.savefig(path + sample + '/venn'+ pontype+ '.png')   # save the figure to file
plt.close(fig) 



### GENERATING VAF & RD HISTOGRAMS
def roundup(x):
    """
    This function rounds up a number to the first higher 200-fold

    Argument: Float or Int

    Returns 200-fold int
    """
    return int(math.ceil(x / 200.0)) * 200

hist_AF = df['AF_mean'].plot_bokeh(
    kind="hist",
    bins=np.linspace(0, 1, 101),
    histogram_type="sidebyside",
    vertical_xlabel=True,
    hovertool=True,
    title="Allele Frequency " + sample,
    ylabel = '# Of mutations',
    xlabel = 'Mean Allele Frequency per mutation, calculated by Mutect2, LoFreq & VarDict',
    legend = None,
    line_color="black",
    show_figure = False)

hist_AF2 = df['AF_mean'].plot_bokeh(
    kind="hist",
    color='darkblue',
    bins=np.linspace(0, 0.01, 21),
    histogram_type="sidebyside",
    vertical_xlabel=True,
    hovertool=True,
    title="Allele Frequency " + sample,
    ylabel = '# Of mutations',
    xlabel = 'Mean Allele Frequency per mutation, calculated by Mutect2, LoFreq & VarDict, Zoomed in on 0 - 0.001',
    legend = None,
    line_color="black",
    show_figure = False)

hist_DP = df['DP_mean'].plot_bokeh(
    kind="hist",
    color='green',
    bins=np.linspace(0, roundup(df['DP_mean'].max()), int(roundup(df['DP_mean'].max())/200+1)),
    histogram_type="sidebyside",
    vertical_xlabel=True,
    hovertool=True,
    title='Read Depth '+ sample,
    ylabel = '# Of mutations',
    xlabel = 'Mean Read Depth per mutation, calculated by Mutect2, LoFreq & VarDict',
    legend = None,
    line_color="black",
    show_figure = False)
    
pandas_bokeh.output_file(path + sample +"/Mean_Depth_&_AF_per_Mutation_"+ pontype + ".html")  
pandas_bokeh.plot_grid([[hist_AF], [hist_AF2], [hist_DP]], plot_width=1500, plot_height=440)


"""Dit_conductance_method.py: Calculating Dit using conductance method."""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import sys

def plotCV(raw_df, name):
    pivot_df = raw_df.pivot(index='VBias', columns='Freq', values='C')
    pivot_df.plot()
    plt.title(name + '_C-V')
    plt.xlabel('V(V)')
    plt.ylabel('C(F)')
    plt.subplots_adjust(right=0.8)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Freq(Hz)', fontsize = 7)
    plt.savefig(name+'_C-V')

def plotGV(raw_df, name):
    pivot_df = raw_df.pivot(index='VBias', columns='Freq', values='G')
    pivot_df.plot()
    plt.title(name + '_G-V')
    plt.xlabel('V(V)')
    plt.ylabel('G(S)')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.subplots_adjust(right=0.8)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Freq(Hz)', fontsize = 7)
    plt.savefig(name+'_G-V')

def plotDit_w(raw_df,r,name):
    df_Cm = raw_df.pivot(index='Freq', columns='VBias', values='C')
    df_Gm = raw_df.pivot(index='Freq', columns='VBias', values='G')
    cox = max(df_Cm.loc[:,max(df_Cm.columns)])  # acc. capacitance of n-type substrate, get the largest acc. capacitance value among all frequencies.
    #cox = max(df_Cm.loc[:,min(df_Cm.columns)])  # acc. capacitance of p-type substrate, get the largest acc. capacitance value among all frequencies.
    df_cox = pd.DataFrame(cox, index=np.arange(len(df_Cm.index)), columns=np.arange(len(df_Cm.columns)))
    df_w = pd.DataFrame(df_Cm.index*2*math.pi)
    for i in range (1,len(df_Cm.columns)):
        df_w[i]=df_Cm.index*2*math.pi
        
    #===== Calculating Gp/w =====#
    data = df_w.values*df_Gm.values*(df_cox**2)/(df_Gm.values**2+(df_w.values**2)*((df_cox-df_Cm.values)**2))
    Gp_w = pd.DataFrame(data=data.values, columns=df_Cm.columns, index=df_Cm.index*2*math.pi)

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(Gp_w.index,Gp_w.values)
    plt.title(name + '_Gp/w-w')
    plt.xscale('log')
    plt.xlabel('w(Hz)')
    plt.ylabel('Gp/w')
    plt.savefig(name + '_Gp_w-w')

    # Plot Gp-w - f for every voltage bias
    # for i in range(0, len(df_Cm.columns)):
     # plt.figure()
     # df.iloc[:,i].plot()

    #===== Calculating Dit =====#
    q = 1.602e-19  # Coulombs of an elementary charge
    Dit = Gp_w.max(axis=0)*2.5/(math.pi*(r**2)*q)
    # Output Dit-V to csv
    Dit.to_csv(name + '_Dit-V', sep='\t', header = ['Dit'])
    print ('Dit = %e (eV-1cm-2)'  %Dit.min(axis=0))
    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)
    ax3.plot(Dit.index,Dit.values)
    plt.title(name + '_Dit-V')
    plt.xlabel('V(V)')
    plt.ylabel('Dit(eV-1cm-2)')
    plt.savefig(name + '_Dit-V')

    return cox

def main():
    raw_df = pd.read_csv(sys.argv[0], sep='\t')  # sys.argv[0] is your your frequency sweep Cm-Gm file path
    name = input('Enter the sample number:')
    r = float(input('Enter the radius of the MOS device(in um):'))*1e-4  #change unit to cm
    d = float(input('Enter the total thickness of the MOS device(in nm):'))*1e-7  #change unit to cm
    print ('Calculating Dit, k and generating figures.....')
    plotCV(raw_df, name)
    plotGV(raw_df, name)
    cox = plotDit_w(raw_df,r, name)
    print ('cox = %e (F)' %cox)
    k0 = 8.85e-14 # vacuum permitivity (F/cm)
    k = float(cox*d/((math.pi*(r**2))*k0))
    print('k = %f' %k)
    plt.show()

if __name__ == '__main__':
    main()

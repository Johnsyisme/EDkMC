import numpy as np
import re
from matplotlib import rcParams
params={'font.family':'sans-serif','font.sans-serif':'Arial','font.style':'normal','font.weight':'normal'}
def extractFile():
    file = "energy_MC.log"
    f2 = open('accept.log','w')
    with open(file, 'r') as f:
        lines = f.readlines()
        for i in lines:
            n = re.findall(r"Accept",i)
            if n:
                print(i)
                f2.writelines(i)
    f2.close()

def monteCarloPlot():
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import MultipleLocator
    import matplotlib.ticker as ticker

    with open("accept.log", 'r') as f:
        lines = f.readlines()
    step = []
    energy = []
    for line in lines:
        tem = re.split(r'[:,]', line)
        step.append(tem[1])
        energy.append(tem[5])
    x=np.array(step,dtype=float)/1000    
    y=np.array(energy,dtype=float)
    n_atom = 744
    y /= n_atom
    plt.figure(dpi=120,figsize=(8,5))
    ax = plt.gca()
    plt.grid(which='major', axis='y',color='lightgray',linestyle='-', linewidth=0.4)
    plt.tick_params(top=False,bottom=True,left=True,right=False,\
        which='major',direction='in',width=3, length=11) 
    Xmax = 200
    plt.scatter(x,y,marker='o',s=1.2,c='crimson')
    plt.hlines(-8.802, 0, Xmax, linestyle='dashed',color='dodgerblue',linewidth=3)
#    plt.annotate("Gr", (12, -9.225),fontsize=24)
    plt.annotate("ma-BN", (5, -5.7), fontsize=24)
    plt.xlim(left=-10, right=Xmax+3)
    plt.xlabel("Step ($\\times$$10^{3}$)",fontsize=24,color="k")
    # plt.xlabel(r'$\mathrm{\mathsf{10^3}}$',fontsize=28,labelpad=6,color="k")
    plt.ylabel('Energy/atom (eV)',labelpad=12,fontsize=24,color="k")
    plt.xticks(np.arange(0,Xmax+1,Xmax/4), fontsize=24)
    plt.yticks(fontsize=24)
    plt.ylim(top = -7.6,bottom = -9.0)
    ymajor = MultipleLocator(0.4)
    ax.yaxis.set_major_locator(ymajor)  
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))  
    #plt.title("45 $\AA$*45 $\AA$ ",fontsize=24)
    # plt.rcParams.update({"text.usetex": True,'font.family':'sans-serif','font.sans-serif':['Arial']})
    plt.savefig("fine_kMC_energy_test.png",dpi=180,transparent=False,bbox_inches='tight')
    plt.tight_layout()
    plt.show()  

extractFile()
monteCarloPlot()



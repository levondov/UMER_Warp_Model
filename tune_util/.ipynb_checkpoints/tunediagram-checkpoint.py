import matplotlib.pyplot as plt
import numpy as np


def tunediagram(order=range(1,4),integer=[0,0],lines=[1,1,1,1],colors='ordered',linestyle='-',fig=plt.gcf()):
    '''
    plot resonance diagram up to specified order
    mx + ny = p
    x = (p-ny)/m
    x = 1 where y = (p-m)/n

    EXAMPLE:
       tunediagram(order=[1,2,3])
       tunediagram([1,2,3],integer=[6,8],lines=[0,0,1,1])
       tunediagram([1,2,3],integer=[6,8],lines=[0,0,1,1], colors='black3', linestyle='--')

    INPUT:
      1. order - array of tune orders to plot. e.g. up to 3rd order, [1,2,3]
      2. integers - integer part of the tune to plot in x,y, default is [0,0]. e.g. plot from 6-7 and 9-10, integer=[6,9]
      2. lines - a boolean of which resonance lines to plot. e.g. [vertical,horizontal,sum,diff]. e.g. plot only vert/horz lines, lines = [1,1,0,0]
      4. colors - option to plot lines in different colors. default is ordered. color options only go up to 10th order
        ordered - each order resonance is a different color
        black - all lines are black
        blackX - X here is a number. all resonances X+ will be in black. e.g. black3 means plot resonances 1-2 in color and 3,4,5,... in black
      6. linestyle - linestyle option from matplotlib
      7. fig - pass in a handle to a figure, otherwise things will be plotted in the current figure.

    Written by Levon Dovlatyan
    University of Maryland, Department of Physics
    Oct 2018
    '''
    # define some variables
    pval = 40 # increase for more/higher order lines
    p = np.linspace(0,pval,pval+1)
    
    qxmin,qymin = integer[0],integer[1]
    qxmax,qymax = qxmin+1,qymin+1
    
    # define different colors, up to 10th order
    color = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9']
    if colors == 'black':
        color = ['k']*10
    elif colors[0:-1] == 'black':
        idx = int(colors[-1])
        color = color[0:idx-1] + (['k']*10)[idx-1:]
    
    # adjust plot limits
    plt.xlim((qxmin-0.01, qxmax+0.01))
    plt.ylim((qymin-0.01, qymax+0.01))
    
    # Plotting formula
    # we plot resonances in reverse order
    for i in order[::-1]:
        m = np.linspace(-i,i,2*i+1)
        n1 = (i-np.abs(m))
        n2 = -1*n1
        for j in range(0,m.size,1):
            # check to see equation is divided by 0 or not
            # ver & hor res lines
            if ((n1[j] == 0 and lines[1]) or (m[j] == 0 and lines[0])):
                # vertical lines
                if n1[j] == 0 and lines[1]:
                    plt.vlines(p/m[j],qymin,qymax,color=color[i-1],linestyle=linestyle);
                # horizontal lines
                if m[j] == 0 and lines[0]:
                    plt.hlines(p/n1[j],qxmin,qxmax,color=color[i-1],linestyle=linestyle);
                    plt.hlines(p/n2[j],qxmin,qxmax,color=color[i-1],linestyle=linestyle);
            # sum and dif res lines
            elif not(n1[j] == 0) and not(m[j] == 0):
                # resonance sum lines
                if lines[2]:
                    if np.sign(m[j]) > 0:
                        plt.plot([[qxmin]*p.size,[qxmax]*p.size],[p/n2[j] - np.array(m[j]*qxmin/n2[j]), p/n2[j] - np.array(m[j]*qxmax/n2[j])],color=color[i-1],linestyle=linestyle);
                    else:
                        plt.plot([[qxmin]*p.size,[qxmax]*p.size],[p/n1[j] - np.array(m[j]*qxmin/n1[j]), p/n1[j] - np.array(m[j]*qxmax/n1[j])],color=color[i-1],linestyle=linestyle);
                # resonance dif lines
                if lines[3]:
                    if np.sign(m[j]) > 0:
                        plt.plot([[qxmin]*p.size,[qxmax]*p.size],[p/n1[j] - np.array(m[j]*qxmin/n1[j]), p/n1[j] - np.array(m[j]*qxmax/n1[j])],color=color[i-1],linestyle=linestyle);
                    else:
                        plt.plot([[qxmin]*p.size,[qxmax]*p.size],[p/n2[j] - np.array(m[j]*qxmin/n2[j]), p/n2[j] - np.array(m[j]*qxmax/n2[j])],color=color[i-1],linestyle=linestyle);

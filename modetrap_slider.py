#!/usr/bin/env python

import matplotlib
matplotlib.use('TkAgg')
from pylab import *
from matplotlib.widgets import Slider, Button, RadioButtons
import modetrap_sub


# Number of beads
nbeads = 1

narg=len(sys.argv)
if narg >=2 :
    nbeads = int(sys.argv[1])


# This is a wrapper for modetrap that insures that the bead locatations ('xpert')
# are sorted into increasing order.
def mode_wrap(n1,n2,xpert,amp,width):
    xp=np.array(xpert)
    ap=np.array(amp)
    wp=np.array(width)
    iargs = np.argsort(xpert)
    xp=xp[iargs]
    ap=ap[iargs]
    wp=wp[iargs]
    res = modetrap_sub.modetrap(n1,n2,xp,ap,wp)
    return res

# Calculate the forward period difference from a set of consecutive periods.
def dp_calc(periods):
    i0 = np.arange(len(periods)-1)
    per = periods[i0]
    dp = periods[i0+1]-periods[i0]
    #return per,dp
    return i0+1,dp

def test_data(i):
    if nbeads==1:
        plists = [ 
            [ 0.27, 0.78, 0.05 ], 
            [ 0.13, 1.44, 0.02 ], 
            [ 0.36, 0.55, 0.03 ],
            [ 0.50, 1.35, 0.04 ],
            [ 0.12, 1.73, 0.07 ] 
            ]
    elif nbeads==2:
        plists = [ 
            [ 0.27, 0.78, 0.05, 0.23, 0.94, 0.03 ], 
            [ 0.10, 1.41, 0.12, 0.13, 1.44, 0.02 ], 
            [ 0.29, 1.32, 0.04, 0.35, 1.50, 0.06 ],
            ]
    else:
        print 'Case for',nbeads,'not supported'
        exit()
    params = plists[i]
    if nbeads==1:
        periods = mode_wrap(n1,n2,[params[0]],[params[1]],[params[2]])
    else:
        periods = mode_wrap(n1,n2,[params[0],params[3]],[params[1],params[4]],[params[2],params[5]])
    per,dp = dp_calc(periods)
    return per,dp

# Overtone/harmonic values "n" for the first ("n1") and last ("n2") mode to be computed.
n1=1
n2=21
# Default/initial parameters for the location, amplitude (mass), and width of the beads.
locvec0 = [ 0.16, 0.22, 0.33 ]
ampvec0 = [ 0.0, 0.0, 0.0 ]
widvec0 = [ 0.03, 0.03, 0.03 ]

# Set the initial set of parameters.
loc = []
amp = []
width = []
for i in np.arange(nbeads):
    loc.append(locvec0[i])
    amp.append(ampvec0[i])
    width.append(widvec0[i])

# Set the initial x and y limits on the dP vs. P plot
#x1lim = np.pi*(n1-1)
#x2lim = np.pi*(n2)
x1lim = n1-0.5
x2lim = n2-0.5
y1lim = 2.0
y2lim = 4.0

xleft = 0.25

fig = figure(1,figsize=(12,6))
ax = subplot(111)
subplots_adjust(left=xleft, bottom=0.35, right=0.975, top=0.95)
periods =  mode_wrap(n1,n2,loc,amp,width)
per,dp = dp_calc(periods)

l, = plot(per,dp, 'ro-', lw=2, color='red')

per0, dp0 = test_data(0)
l0, = plot(per0,dp0, 'ko--', lw=3, color='black', markersize=10, mfc='none')

axis([x1lim, x2lim, y1lim, y2lim])
#xlabel(r'$\omega_n\, \rm (frequency)$', fontsize=22)
xlabel(r'$n \,\, \rm (overtone\, number)$', fontsize=22)
ylabel(r'$\omega_{n+1} \,-\, \omega_n$',fontsize=22)

axcolor = 'lightgoldenrodyellow'

axvec = []
slvec = []
savec = []
swvec = []
xright = 0.93
xspace = 0.1
h = ((xright-xleft) - (nbeads-1)*xspace)/real(nbeads)
for i in np.arange(nbeads):
    xa = xleft + real(i)*(h+xspace)
    lab = 'Bead {0}'.format(i+1)
    axamp  = axes([xa, 0.15, h, 0.03], axisbg=axcolor)
    text(1.2,1.50,lab)
    axloc  = axes([xa, 0.10, h, 0.03], axisbg=axcolor)
    axwid  = axes([xa, 0.05, h, 0.03], axisbg=axcolor)
    axs = [ axloc, axamp, axwid ]
    axvec.append(axs)
    if i==0:
        sloc = Slider(axloc, 'Location ',0.0, 1.0, valinit=locvec0[i])
        samp = Slider(axamp, 'Mass ',    0.0, 3.0, valinit=ampvec0[i])
        swid = Slider(axwid, 'Width ',   0.0, 0.4, valinit=widvec0[i])
    else:
        sloc = Slider(axloc, '', 0.0, 1.0, valinit=locvec0[i])
        samp = Slider(axamp, '', 0.0, 3.0, valinit=ampvec0[i])
        swid = Slider(axwid, '', 0.0, 0.4, valinit=widvec0[i])
    slvec.append(sloc)
    savec.append(samp)
    swvec.append(swid)

def update(val):
    ampvec = []
    locvec = []
    widvec = []
    for i in np.arange(nbeads):
        ampvec.append(savec[i].val)
        locvec.append(slvec[i].val)
        widvec.append(swvec[i].val)
    periods =  mode_wrap(n1,n2,locvec,ampvec,widvec)
    per,dp = dp_calc(periods)
    l.set_xdata(per)
    l.set_ydata(dp)
    draw()

for i in np.arange(nbeads):
    slvec[i].on_changed(update)
    savec[i].on_changed(update)
    swvec[i].on_changed(update)

resetax = axes([0.025, 0.100, 0.075, 0.06])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
def reset(event):
    for i in np.arange(nbeads):
        slvec[i].reset()
        savec[i].reset()
        swvec[i].reset()
button.on_clicked(reset)

rax = axes([0.025, 0.3, 0.12, 0.50], axisbg=axcolor, aspect='equal')
blabs = []
for i in np.arange(1,ncases+1):
    blabs.append('Case {0}'.format(i))

#radio = RadioButtons(rax, ['Case 1', 'Case 2', 'Case 3'], active=0)
radio = RadioButtons(rax, blabs, active=0)
def colorfunc(label):
    datakey = int(label[-2:])
    per0, dp0 = test_data(datakey-1)
    l0.set_xdata(per0)
    l0.set_ydata(dp0)
    #l.set_color(label)
    draw()
radio.on_clicked(colorfunc)

show()


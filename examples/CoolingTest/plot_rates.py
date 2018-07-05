import matplotlib
matplotlib.use("Agg")
from pylab import *
from scipy import stats
import h5py

# Plot parameters
params = {'axes.labelsize': 10,
'axes.titlesize': 10,
'font.size': 12,
'xtick.labelsize': 10,
'ytick.labelsize': 10,
'text.usetex': True,
 'figure.figsize' : (3.15,2.65),
'figure.subplot.left'    : 0.19,
'figure.subplot.right'   : 0.99,
'figure.subplot.bottom'  : 0.14,
'figure.subplot.top'     : 0.99,
'figure.subplot.wspace'  : 0.15,
'figure.subplot.hspace'  : 0.12,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})

# Read the file header 
file = open("rates.txt", "r")
line_redshift = file.readline()
line_abundance = file.readline()
file.close()

# This is the most dirty line of code
redshift = float(line_redshift[13:])
abundance = float(line_abundance[12:])

# Read the rates for the different densities 
data = loadtxt("rates.txt")
T_m6 = data[:,0]
du_dt_m6 = data[:,1]
T_m5 = data[:,2]
du_dt_m5 = data[:,3]
T_m4 = data[:,4]
du_dt_m4 = data[:,5]
T_m3 = data[:,6]
du_dt_m3 = data[:,7]
T_m2 = data[:,8]
du_dt_m2 = data[:,9]
T_m1 = data[:,10]
du_dt_m1 = data[:,11]
T_0 = data[:,12]
du_dt_0 = data[:,13]

#Plot limis
x_min = 3e3
x_max = 5e8
y_min = 1e-24
y_max = 3e-21

# Plot the net heating and cooling rates for a selection of densities
figure()
loglog(T_m6, du_dt_m6, '--', color='C0', lw=1)
loglog(T_m6, -du_dt_m6, '-', color='C0', lw=1, label="$n_{\\rm H} = 10^{-6}$")

loglog(T_m4, du_dt_m4, '--', color='C1', lw=1)
loglog(T_m4, -du_dt_m4, '-', color='C1', lw=1, label="$n_{\\rm H} = 10^{-4}$")

loglog(T_m2, du_dt_m2, '--', color='C2', lw=1)
loglog(T_m2, -du_dt_m2, '-', color='C2', lw=1, label="$n_{\\rm H} = 10^{-2}$")

loglog(T_0, du_dt_0, '--', color='C3', lw=1)
loglog(T_0, -du_dt_0, '-', color='C3', lw=1, label="$n_{\\rm H} = 10^{0}$")

text(3.5e8, 2.35e-21, "$z=%.3f$\n$Z/Z_\\odot=%.3f$"%(redshift, abundance), ha="right", va="top", fontsize=8.5, backgroundcolor='w')

legend1 = legend(loc="lower right", ncol=2, frameon=False,  handletextpad=0.2, handlelength=0.8, fontsize=8.5)

plot_lines = []
l1, = plot([-1, -1], 'k--', lw=1.)
l2, = plot([-1, -1], 'k-', lw=1.)
plot_lines.append(l1)
plot_lines.append(l2)
legend2 = legend(plot_lines, ["${\\rm Heating}$", "${\\rm Cooling}$"], loc="upper left", frameon=False,  handletextpad=0.2, handlelength=1.2, fontsize=8.5)

xlim(x_min, x_max)
ylim(y_min, y_max)
xlabel("$T~[{\\rm K}]$", labelpad=0)
ylabel("$|\Lambda / n_{\\rm H}^2|~[{\\rm erg\\cdot cm^{3}\\cdot s^{-1}}]$", labelpad=-1)
gca().add_artist(legend1)
gca().add_artist(legend2)

savefig("cooling_rates.pdf")

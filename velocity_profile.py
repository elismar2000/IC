from astropy.io import fits
import matplotlib.pyplot as plt

file1 = 'sft00.fits'

#hdul = fits.open(file1)

with fits.open(file1) as hdul:
#    hdul.info()
    hdr = hdul[0].header
    merger = hdul[1].data

x = []
v_y = []
for i in range(merger.shape[0]):
    if abs(merger[i][2]) < 0.005 and merger[i][0] > 50.0 and merger[i][0] < 80.0:
        x.append(merger[i][0])
        v_y.append(merger[i][4])

plt.plot(x,v_y,'b.',label='')
plt.xlabel('x (kpc)')
plt.ylabel('v_y (m/s)')
plt.show()

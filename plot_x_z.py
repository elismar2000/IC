from astropy.io import fits
import matplotlib.pyplot as plt

file1 = 'sft00.fits'

#hdul = fits.open(file1)

with fits.open(file1) as hdul:
#    hdul.info()
    hdr = hdul[0].header
    merger = hdul[1].data

x = []
z = []
for i in range(merger.shape[0]):
    x.append(merger[i][0])
    z.append(merger[i][2])

plt.plot(x,z,'b.')
plt.show()

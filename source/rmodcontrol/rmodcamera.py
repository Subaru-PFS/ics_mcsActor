import subprocess as sub
p = sub.Popen(['/opt/EDTpdv/initcam', '-f', '/home/chyan/71MP/illunis-71mp.cfg'],stdout=sub.PIPE,stderr=sub.PIPE)
output, errors = p.communicate()
print errors[23:-1]

expTime =1000
p = sub.Popen(['rmodcontrol', '-e', str(expTime)],stdout=sub.PIPE,stderr=sub.PIPE)
output, errors = p.communicate()
print output[28:-6]

f = pyfits.open('/home/chyan/mhs/data/mcs/schmidt_fiber_snr400_rmod71.fits')   
#!/usr/bin/python

import struct

prop_main  = 2650.0
prop_layer = 2160.0

# size of the domain
sx = 1200.0
sy = 1200.0

# coordinates of the layer
Y0 = 600.0
Y1 = 1200.0

# grid
nx = 400
ny = 400

tol = 1e-8

hx = sx / nx; print 'hx = ', hx
hy = sy / ny; print 'hy = ', hy


prop_values = []
for iy in range(ny):
  y = (iy+0.5)*hy
  for ix in range(nx):
    if y > Y0 and y < Y1:
      prop_values.append(prop_layer)
    else:
      prop_values.append(prop_main)


# write the properties into a binary file with a single precision
out = open('properties.bin', 'wb')
s = struct.pack('f'*len(prop_values), *prop_values)
out.write(s)
out.close()


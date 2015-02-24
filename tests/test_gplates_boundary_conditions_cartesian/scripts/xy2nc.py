#!/usr/bin/env python

# Script that converts GPLATES output (GMT format) into NetCDF format, compatible with Fluidity.

from Scientific.IO.NetCDF import NetCDFFile
import os.path
import sys
import numpy

def main(input_file_name):
  root, ext = os.path.splitext(input_file_name)
  data      = numpy.loadtxt(input_file_name)
  lon       = numpy.unique(data[:,0])
  lat       = numpy.unique(data[:,1])
  dlon      = lon[1]-lon[0]
  m         = int(round(360./dlon))
  dlat      = lat[1]-lat[0]
  n         = int(round(180./dlat))

  # Some asserts to verify that gplates data is of the dimension expected:
  assert( m == len(lon) )
  assert( n+1 == len(lat) )
  assert( m*(n+1)== data.shape[0] )
  assert( lon[-1]== 179 )

  data = data.reshape( (n+1, m, 6) )
  assert( all(data[1,:,0]    == numpy.linspace(-180.,180., m+1)[:-1]) )
  assert( all(data[::-1,0,1] == numpy.linspace(-90,90., n+1)) )

  nc = NetCDFFile(root+'.nc', 'w')
  nc.createDimension('lon', len(lon)+1)
  nc.createDimension('lat', len(lat))
  nc.createDimension('dim', 3)

  lonvar = nc.createVariable('lon', data.dtype.char, ('lon',))
  lonvar[:] = numpy.linspace(-180., 180., m+1)
  lonvar.long_name = "longitude"
  lonvar.units = "degrees_east"
  lonvar.actual_range = (-180., 180.);

  latvar = nc.createVariable('lat', data.dtype.char, ('lat',))
  latvar[:] = numpy.linspace(-90., 90., n+1)
  latvar.long_name = "latitude"
  latvar.units = "degrees_north"
  latvar.actual_range = (-90, 90.);

  velvar = nc.createVariable('velocity', data.dtype.char, ('dim', 'lon', 'lat'))
  velvar[:,:m,:] = data[::-1,:,3:].transpose((2,1,0)) / (365*24*60*60*100.)
  velvar[:,m,:] = data[::-1,0,3:].transpose() / (365*24*60*60*100.)
  velvar.long_name = "velocity"
  velvar.units = "cm/yr"

  plateid = nc.createVariable('plate_id', data.dtype.char, ('lon', 'lat'))
  plateid[:m,:] = data[::-1,:,2].transpose((1,0))
  plateid[m,:] = data[::-1,0,2]
  plateid.long_name = "plate_id"

  # global attributes copied from gmt tools output
  nc.Conventions = "COARDS/CF-1.0"
  nc.title  = root+'.nc'
  nc.history = " ".join(sys.argv)
  nc.close()

def usage():
  sys.stderr.write('USAGE: xy2nc.py <input_xy_file>\n')

if __name__ == "__main__":
  if len(sys.argv)!=2:
    usage()
    sys.exit(1)
  main(sys.argv[1])

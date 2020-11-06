import gamma.friction_coupling as fc
import sys

a = fc.friction_gamma_parser(aims_file=sys.argv[1],gamma_files=sys.argv[2:])

b = fc.friction_tensor(a.ks,a.coords,a.eis,a.ejs,a.couplings,a.kweights,a.chem_pot,a.friction_masses,300,0.6,nspin=1)

tensor = b.calc_tensor()

print(tensor)
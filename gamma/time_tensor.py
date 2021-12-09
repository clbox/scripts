import gamma.friction_coupling as fc
import glob


dir_name1 = '1/'
dir_name2 = '2/'


gamma_files1 = glob.glob(dir_name1+'/*friction_gamma*.out')
gamma_files2 = glob.glob(dir_name2+'/*friction_gamma*.out')

parser = fc.friction_output_parser_2021()
eigenvalues1, couplings1, n_excits1, k_points1 = parser.parse_gamma_couplings(gamma_files1)
eigenvalues2, couplings2, n_excits2, k_points2 = parser.parse_gamma_couplings(gamma_files2)

chem_pot1 = parser.parse_chem_pot(dir_name1+'aims.out')
chem_pot2 = parser.parse_chem_pot(dir_name2+'aims.out')
print('Chemical potentials / eV : '+ str(chem_pot1) + ' ' + str(chem_pot2) + ' Difference / eV : ' +str(chem_pot1-chem_pot2))
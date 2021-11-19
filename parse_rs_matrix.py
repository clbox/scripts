import numpy as np


class rs_parser():

    def __init__(self):
        print('start')


    def parse_indices(self,indices_file):

        self.column_index_hamiltonian = []

        read_cell_index = False
        read_index_hamiltonian1 = False
        read_index_hamiltonian2 = False
        read_column_index = False


        with open(indices_file,'r') as f:
            for line in f:
                if 'n_hamiltonian_matrix_size:' in line:
                    self.n_hamiltonian_matrix_size = int(line.split()[-1])
                if 'n_cells_in_hamiltonian:' in line:
                    self.n_cells_in_hamiltonian = int(line.split()[-1])
                if 'n_basis:' in line:
                    self.n_basis = int(line.split()[-1])


                if 'column_index_hamiltonian' in line:
                    read_index_hamiltonian2 = False
                    read_column_index = True
                    continue

                if read_column_index:
                   self.column_index_hamiltonian.append(int(line))


                # Begin parse index_hamiltonian1
                if 'index_hamiltonian(2,:,:)' in line:
                    read_index_hamiltonian1 = False
                    read_index_hamiltonian2 = True
                    i_cell = 0
                    continue

                if read_index_hamiltonian2:
                    self.index_hamiltonian[1,i_cell,:] = np.fromstring(line,sep=' ')
                    i_cell +=1
                    continue


                # Begin parse index_hamiltonian1
                if 'index_hamiltonian(1,:,:)' in line:
                    read_cell_index = False
                    read_index_hamiltonian1 = True
                    self.index_hamiltonian = np.zeros((2,self.n_cells_in_hamiltonian,self.n_basis),dtype=np.int64)
                    i_cell = 0
                    continue

                
                if read_index_hamiltonian1:
                    self.index_hamiltonian[0,i_cell,:] = np.fromstring(line,sep=' ')
                    i_cell +=1
                    continue

                if 'cell_index' in line:
                    read_cell_index = True
                    self.cell_index = np.zeros((self.n_cells_in_hamiltonian,3),dtype=np.int64)
                    i_cell = 0
                    continue

                # Begin parse cell index
                if read_cell_index:
                    self.cell_index[i_cell,0] = int(line.split()[0])
                    self.cell_index[i_cell,1] = int(line.split()[1])
                    self.cell_index[i_cell,2] = int(line.split()[2])
                    i_cell += 1

            self.column_index_hamiltonian = np.array(self.column_index_hamiltonian,dtype=np.int64)

    def parse_rs_matrix(self,rs_matrix_file):

        rs_matrix = np.loadtxt(rs_matrix_file)

        return rs_matrix


    def construct_dense_from_sparse(self,rs_matrix,k_phase,i_k_point):

        matrix_complex = np.zeros((self.n_basis,self.n_basis),dtype=np.complex128)
        for i_cell in range(1,self.n_cells_in_hamiltonian):
            for i_basis_1 in range(1, self.n_basis+1):
                if( self.index_hamiltonian[0,i_cell-1, i_basis_1-1] > 0 ): 
                    for i_place in range(self.index_hamiltonian[0,i_cell-1, i_basis_1-1], \
                        self.index_hamiltonian[1,i_cell-1, i_basis_1-1]):

                        i_basis_2 =  self.column_index_hamiltonian[i_place-1]

                        matrix_complex[i_basis_1-1, i_basis_2-1] =  \
                        matrix_complex[i_basis_1-1, i_basis_2-1] +  \
                        np.conj(k_phase[i_cell-1,i_k_point-1])* \
                        rs_matrix[i_place-1]


                        if(i_basis_1!=i_basis_2):
                            matrix_complex[i_basis_2-1, i_basis_1-1] =  \
                            matrix_complex[i_basis_2-1, i_basis_1-1] +  \
                            k_phase[i_cell-1,i_k_point-1]* \
                            rs_matrix[i_place-1]

        return matrix_complex




a = rs_parser()

a.parse_indices('rs_indices.out')

ham = a.parse_rs_matrix('rs_hamiltonian.out')

# i_cell, i_k_point
k_phase = np.zeros((2,8),dtype=np.complex128)
k_phase[0,:] = complex(1.0,0.0)
k_phase[1,:] = complex(0.0,0.0)

print(a.construct_dense_from_sparse(ham,k_phase,1))
                    







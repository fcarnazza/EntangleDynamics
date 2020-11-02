from subroutines import *
run           = 'statevector' # 'qasm' 'ibmq_...'
x0            = [1,1,1,0]
out_RO        = open('1110','w')
coupling      = 0.1
L             = Ising_Lattice(n=4,bonds=[(0, 1), (1, 2), (2,3)])
H             = Ising_Hamiltonian(lattice=L,B=coupling)#
algorithm  = Quench([0,2],H,t0=0.,t=1,tstep=0.1,x0=x0)
algorithm.generate_and_run(run,print_info=True,shots=8192,layout=[0,1,2,3],optimization=3,calibration=True,output=out_RO)

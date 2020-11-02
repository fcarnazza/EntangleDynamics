from qiskit.providers.aer.noise                        import NoiseModel
from qiskit import *
import numpy as np
from scipy import linalg as LA
from qiskit.quantum_info.synthesis.two_qubit_decompose import TwoQubitBasisDecomposer
from qiskit.extensions.standard.x import CnotGate
from qiskit.aqua import QuantumInstance,aqua_globals
from qiskit.aqua.operators import WeightedPauliOperator
from qiskit.quantum_info import Pauli, Operator
from qiskit.ignis.mitigation.measurement import (complete_meas_cal,tensored_meas_cal,CompleteMeasFitter,TensoredMeasFitter)
import itertools,functools
import more_itertools as mit
from math import cos, sin, sqrt, pi, acos, exp, pi
import string

def generate_pauli_Z(idx,n):
    zeros = [0]*n
    zmask = [0]*n
    for i in idx: zeros[i] = 1
    a_x = np.asarray(zeros,dtype=np.bool)
    a_z = np.asarray(zmask,dtype=np.bool)
    return Pauli(a_x,a_z)

def generate_number_op(x,idx,n):
        l = len(idx)
        llist = list(range(l))
        N_list = []
        for i in range(l+1):
                         comb = list(itertools.combinations(idx,i))
                         n_comb = list(map(lambda elem: ((-1)**len([*elem])*x/2**(l),generate_pauli_Z([*elem],n)),comb))
                         N_list.append( n_comb)
        #flattened_N_list = list(itertools.chain(*N_list))
        return N_list

def print_unitary(c,prin=False):
    simulator = Aer.get_backend('unitary_simulator')
    result    = execute(c,simulator).result()
    unitary   = result.get_unitary(c)
    U = np.zeros((16,16), dtype=complex) 
    if(prin):    
            for i in range(unitary.shape[0]):
                for j in range(unitary.shape[1]):
                    if(np.abs(unitary[i,j])>1e-6): print(i,j,unitary[i,j])
    return unitary
    
def transpilation(c,machine,opt,layout):
    provider = IBMQ.load_account()
    provider = IBMQ.get_provider(hub='ibm-q')
    backend  = provider.get_backend(machine)
    from qiskit.compiler import transpile
    transpiled_circuit = transpile(c,backend,initial_layout=layout,optimization_level=opt)
    return transpiled_circuit

def generate_initial_circuit(x0):
    n = len(x0)
    circuit = QuantumCircuit(n,n)
    for i in range(n):
        if(x0[i]==1): circuit.x(i)
    return circuit

def GHZ(n):
    circuit = QuantumCircuit(n,n)
    for control in range(n-1):
      target = control + 1
      circuit.cx(control,target)
    circuit.h(n-1)
    for target in reversed(range(1,n)):
      control = target - 1
      circuit.cx(control,target)
    return circuit

def generate_pauli_dictionary():
    pauli_dict = {}
    pauli_dict['I'] = np.zeros((2,2),dtype=complex)
    pauli_dict['X'] = np.zeros((2,2),dtype=complex)
    pauli_dict['Y'] = np.zeros((2,2),dtype=complex)
    pauli_dict['Z'] = np.zeros((2,2),dtype=complex)
    pauli_dict['I'][0,0] = pauli_dict['I'][1,1] = 1.0
    pauli_dict['X'][0,1] = pauli_dict['X'][1,0] = 1.0
    pauli_dict['Z'][0,0] =  1.0
    pauli_dict['Z'][1,1] = -1.0
    pauli_dict['Y'][0,1] = -1j
    pauli_dict['Y'][1,0] =  1j
    # pauli I
    # 1  0
    # 0  1
    # pauli X
    # 0  1
    # 1  0
    # pauli
    # 0 -i
    # i  0
    # pauli Z
    # 1  0
    # 0 -1
    return pauli_dict

def get_matrices(t,c):
    # 2  0  0  c
    # 0  0  c  0
    # 0  c  0  0
    # c  0  0 -2
    pauli_dict = generate_pauli_dictionary()
    H = np.zeros((4,4))
    H =   np.kron(pauli_dict['Z'],pauli_dict['I']) \
      +   np.kron(pauli_dict['I'],pauli_dict['Z']) \
      + c*np.kron(pauli_dict['X'],pauli_dict['X'])
    U = LA.expm(-1j*t*H)
    #print("Ising matrices ")
    #for r in range(4):
    #    for c in range(4):
    #        print(r,c,H[r,c],U[r,c])
    #print("------")
    return H,U


def tomography(logfile,circuit,quantum_instance):i

    #----------------------------------------------------------------------------------------------------
    # At this point non scalability is apparent: one measures 256 expectation values!
    # This part need to be improved to reach scalability along the lines of e.g. https://arxiv.org/pdf/cond-mat/0503393.pdf
    #----------------------------------------------------------------------------------------------------
    paulis    = [[0,0],[0,1],[1,1],[1,0]]
    names     = ['I','X','Y','Z']
    names_mn  = []
    paulis_mn = []
    #logfile.write("Tomography")
    i,j,k,l = 0,1,2,3

    #provider     = IBMQ.load_account()

    # generiamo gli operatori di Pauli per 2 qubit
    for m,(ax,az) in enumerate(paulis):
        for n,(bx,bz) in enumerate(paulis):
            for p,(cx,cz) in enumerate(paulis):
                for q,(dx,dz) in enumerate(paulis):
                    x_vector = [0]*circuit.n_qubits
                    z_vector = [0]*circuit.n_qubits
                    x_vector[i] = ax
                    z_vector[i] = az
                    x_vector[j] = bx
                    z_vector[j] = bz
                    x_vector[k] = cx
                    z_vector[k] = cz
                    x_vector[l] = dx
                    z_vector[l] = dz
                    x_vector    = np.asarray(x_vector,dtype=np.bool)
                    z_vector    = np.asarray(z_vector,dtype=np.bool)
                    Pmn         = WeightedPauliOperator([(1.0,Pauli(x_vector,z_vector))])
                    paulis_mn.append(Pmn)
                    names_mn.append(names[m]+names[n]+names[p]+names[q])

    #----------------------------------------------------------------------------------------------------
    #Here all circuits to evaluate Pauli operator products are built:
    #----------------------------------------------------------------------------------------------------

    circuits = []
    for idx,oper in enumerate(paulis_mn):
        c = oper.construct_evaluation_circuit(wave_function=circuit,
            statevector_mode=quantum_instance.is_statevector,
            use_simulator_snapshot_mode=False,
            circuit_name_prefix=str(idx))
        circuits.append(c)
     #To print all the details of this circuits, uncomment the following lines:
        #IBMQ.load_account()
        #provider = IBMQ.get_provider(hub='ibm-q',group='open',project='main')
        #backend      = provider.get_backend('ibmqx2')
        #from qiskit.compiler import transpile
        #transpiled_circuit = transpile(c[0],backend,optimization_level=0)
        #print("circuit for operator ",names_mn[idx])
        #print(transpiled_circuit.draw()) 

    #----------------------------------------------------------------------------------------------------
    # Expectations values for this products are computed: 
    #----------------------------------------------------------------------------------------------------

    if circuits:
        to_be_simulated_circuits = \
            functools.reduce(lambda x, y: x + y, [c for c in circuits if c is not None])
        result = quantum_instance.execute(to_be_simulated_circuits)

    #----------------------------------------------------------------------------------------------------
    # And stored on a dictionary: 
    #----------------------------------------------------------------------------------------------------

    tomography_dict = {}
    for idx,oper in enumerate(paulis_mn):
        mean, std = oper.evaluate_with_result(
                    result=result,statevector_mode=quantum_instance.is_statevector,
                    use_simulator_snapshot_mode=False,
                    circuit_name_prefix=str(idx))
        mean,std = np.real(mean),np.abs(std)
        tomography_dict[names_mn[idx]] = [mean,std]
    for k in tomography_dict.keys():
        logfile.write("Pauli operator %s %s\n"%(k, tomography_dict[k]))
    return tomography_dict

def minus_x_log_x(s):
    if(np.abs(s)<1e-6): return 0.0
    return - s*np.log(s)

def analyze_and_print(logfile,rho,name):
    rho /= rho.shape[0]
    #logfile.write("Density operator "+str(name))
    #for r in range(rho.shape[0]):
    #    for c in range(rho.shape[1]):
    #        rh=rho[r,c]
            #logfile.write("rho[%d,%d]=%f \n" % (r,c,rho[r,c]))
    #logfile.write("Trace %f \n" % np.trace(rho))
    #logfile.write("Hermiticity %f \n" % np.abs(rho-np.conj(rho.T)).max())
    sigma,U = LA.eigh(rho)
    #rho_tilde = np.einsum('ik,k,kl->il',U,sigma,np.conj(U.T))

    #----------------------------------------------------------------------------------------------------
    # The density matrix is reconstructed from rho:
    #----------------------------------------------------------------------------------------------------

    #sigma   = np.abs(sigma)
    sigma   = [ max(s,0) for s in sigma]
    sigma  /= np.sum(sigma)
    #logfile.write("Eigenvalues \n"+str(sigma))
    P       = np.sum(sigma**2)
    S       = np.sum([minus_x_log_x(s) for s in sigma])
    #logfile.write("Purity, von Neumann entropy %f, %f \n" %(P,S))

    rho_tilde = np.einsum('ik,k,kl->il',U,sigma,np.conj(U.T))

    return P,S,rho_tilde

# ====================================================================== #

def process_tomography_dict(logfile,tomography_dict,idx=None,two_sigma_truncation=True):

    #----------------------------------------------------------------------------------------------------
    # For a system of k qubits:
    # rho = (\sum_m \langle P_m \rangle P_m)/2**k
    # where  m is the multi-index m_1,...m_n with m_i \in (0,1,2,3)
    # e.g. for just one qubit:
    # rho = (I+a*X+b*Y+c*Z)/2
    # for (a,b,c) \in S^2 the Bloch sphere
    #----------------------------------------------------------------------------------------------------
    
    pauli_dict = generate_pauli_dictionary()
    n_sample = 100
    I_vector = S_A_vector = S_B_vector= S_AB_vector = np.zeros(n_sample)

    #----------------------------------------------------------------------------------------------------
    # To retrieve error bars, for such deep circuits, the best choice is
    # to repeat some times the expertiments. Since such experiments are lengnthy, a
    # possible, but not always satisfying procedure, is to statistically generate them:
    #----------------------------------------------------------------------------------------------------

    for i in range(n_sample):
        rho_total = np.zeros((16,16),dtype=complex)# 4qb->(16,16)
        for k in tomography_dict.keys():
            p0,p1   = k[0],k[1]
            p2,p3  = k[2],k[3]
            ave,std = tomography_dict[k]
            if(two_sigma_truncation):
               if(np.abs(ave)>2*std):
                  sample  = np.random.normal(ave,std,1)
                  krp0p1  = np.kron(pauli_dict[p0]/2,pauli_dict[p1]/2)
                  krp2p3  = np.kron(pauli_dict[p2]/2,pauli_dict[p3]/2)
                  rho_total += sample*np.kron(krp0p1,krp2p3)
            else:
               sample  = np.random.normal(ave,std,1)
               krp0p1  = np.kron(pauli_dict[p0]/2,pauli_dict[p1]/2)
               krp2p3  = np.kron(pauli_dict[p2]/2,pauli_dict[p3]/2)
               rho_total += sample*np.kron(krp0p1,krp2p3)
        n,m = 4,8
        x = [t for t in string.ascii_lowercase[:m]]
        y = [t for t in string.ascii_lowercase[:m]]
        #to_remove = []
        #print(''.join(x))
        for k in range(n):
                if(k not in idx):
                        y.remove(x[k]); y.remove(x[n+k])
                        x[n+k] = x[k]

        # rho(I,J) = <i1 i2 ... iN-1 iN|rho|j1 j2 ... jN-1 jN>
        rho_total = rho_total.reshape((2,2,2,2,2,2,2,2))
        rho_AB = np.einsum(''.join(x)+'->'+''.join(y),rho_total)
        #print("x",x,"y",y)
        rho_AB  = rho_AB.reshape((2,2,2,2))
        rho_AB  = rho_AB.reshape((4,4))
        P_AB,S_AB,rho_AB = analyze_and_print(logfile,rho_AB,'rho(A,B)')
        rho_AB  = rho_AB.reshape((2,2,2,2))
        rho_A  = np.einsum('ikjk->ij',rho_AB)#...
        rho_A  = rho_A.reshape((2,2))
        rho_B  = np.einsum('mimj->ij',rho_AB)#...
        rho_B  = rho_B.reshape((2,2))
        P_A,S_A,_ = analyze_and_print(logfile,rho_A,'rho(A)')
        P_B,S_B,_ = analyze_and_print(logfile,rho_B,'rho(B)')
        I         = S_B+S_A-S_AB
        I_vector[i] = I
        #if(I==0):
        rho_AB  = rho_AB.reshape((4,4))
        #print("rho_AB",rho_AB,"rho_A",rho_A, "rho_B", rho_B )
        #e, u= LA.eigh(rho_AB)
        #print("rho eigenvalues", e)
        #eA, u= LA.eigh(rho_A)
        #print("rho eigenvalues", e)
        #eB, u= LA.eigh(rho_B)
        #print("rhoAB eigenvalues %s, rhoA eigenvalues %s, rhoB eigenvalues %s "%(e,eA,eB) )
    Iave,Istd = np.mean(I_vector),np.std(I_vector)
    S_Aave, S_Astd = np.mean(S_A_vector),np.std(S_A_vector)
    S_Bave, S_Bstd = np.mean(S_B_vector),np.std(S_B_vector)
    S_ABave, S_ABstd = np.mean(S_AB_vector),np.std(S_AB_vector)
    return Iave,Istd,S_Aave,S_Astd,S_Bave,S_Bstd,S_ABave,S_ABstd, rho_AB, rho_A, rho_B


# ===============================================================================================

class Ising_Lattice:
      def __init__(self,n=None,bonds=None):#,logfile=None):
          self.n     = n
          self.bonds = bonds
      def eo_bonds(self):
          even_bonds = []
          odd_bonds = []
          for (i,j) in self.bonds:
              if(i % 2 == 0):
                even_bonds.append((i,j))
              else:
                        odd_bonds.append((i,j))
          return even_bonds, odd_bonds
          #if(logfile is not None):
            # logfile.write("lattice sites %d\n" % self.n)
             #for (i,j) in self.bonds:
            #     logfile.write("dimer: %d %d\n" % (i,j))


# ===============================================================================================

class Ising_Hamiltonian:
        def __init__(self,lattice=None,B=0):#,logfile=None):
                  self.B       = B
                  self.lattice = lattice
#==============================================================================================
# functions: evolve_all_4_spins_up, evolve_dimer, evolve_odd, evolve_even, where employed only during debugging, and may be deleted
#==============================================================================================
        def evolve_all_4_spins_up(self,time):                   #
            Ht = 2*time*sqrt(1+self.B**2)                       #   
            cir = QuantumCircuit(4,4)                           #  here one change the basis from the spin basis state in which  
            phi = -np.sign(B - cos(2*pi*labels[i]/n))*acos(self.B/sqrt(1+self.B**2))   #  all spins are up: |0000>, to the basis of the diagonalized H: 
            cir.rx(phi,2)                                       #  U    |0000> = |00>(cos(phi)|00> + i sin(phi)|11>)  
            #circuit.u1(Ht,0)                                   #   dis   
            cir.cx(2,3)                                         #  and then evolves this state according to the free Hamiltonian.
            cir = self.evolve_diagonal(cir,time)
            return cir                                          #   
        def evolve_dimer(self,qc,i,j,dt):                       #
                c=self.B
                H,U  = get_matrices(dt,c)
                two_qubit_cnot_decompose = TwoQubitBasisDecomposer(CnotGate())
                C = two_qubit_cnot_decompose.__call__(U)
                parameter_string = []
                for g in C:
                    instruction,q1,q2 = g
                    if(instruction.name=='u3'):
                       t1,t2,t3 = instruction.params
                       for x in [t1,t2,t3]: parameter_string.append(round(x,4))
                       if(q1[0].index==0): idx = i
                       else:               idx = j
                       qc.u3(t1,t2,t3,idx)
                    if(instruction.name=='cx'):
                       if(q1[0].index==0): idx_ctrl = i; idx_targ = j
                       else:               idx_ctrl = j; idx_targ = i
                       qc.cx(idx_ctrl,idx_targ)
                return qc
        def evolve_odd(self,circuit,dt,barr=False):
                e,o =self.lattice.eo_bonds()
                for i,j in o:
                    self.evolve_dimer(circuit,i,j,dt)
                if(barr): circuit.barrier()
                return circuit

        def evolve_even(self,circuit,dt,barr=False):
                e,o =self.lattice.eo_bonds()
                for i,j in e:
                    self.evolve_dimer(circuit,i,j,dt)
                if(barr): circuit.barrier()
                return circuit

#==============================================================================================


        def fSWAP(self,qc,i,j):
                qc.swap(i,j)
                qc.cz(i,j)
                return qc
                 
        def Ph_Delay_Beam_Spl(self,qc,qubit_i,qubit_j,k, dag = False):  #    
                n = self.lattice.n                                      #     
                if(dag==True):                                          #        
                        qc.cz(qubit_j,qubit_i)                          #
                        qc.cx(qubit_j,qubit_i)                          #                             n
                        qc.ch(qubit_i,qubit_j)                          #       This circuit is the  F  :  (https://arxiv.org/pdf/1807.07112.pdf)
                        qc.cx(qubit_j,qubit_i)                          #                             k
                        qc.u1(-2*pi*k/n,qubit_j)                        #               
                if(dag==False):                                         #       
                        qc.u1(2*pi*k/n,qubit_j)                         #       |i>----------------------[X]----[X]--[Z]-- 
                        qc.cx(qubit_j,qubit_i)                          #                                 |   |   |    |
                        qc.ch(qubit_i,qubit_j)                          #                                 |   |   |    |
                        qc.cx(qubit_j,qubit_i)                          #       |j>--------U1(2 pi k/n)-----[H]---------      
                        qc.cz(qubit_j,qubit_i)                                 
                return qc                        
                 
        def Bogoliubov(self,qc,qubit_i,qubit_j,k,dag = False):
                B = self.B
                n = self.lattice.n
                theta = acos(-(B - cos(2*pi*k/n)) /sqrt((B - cos(2*pi*k/n))**2 + (sin(2*pi*k/n ))**2 ))
                #print("theta",theta/2)
                qc.x(qubit_i)
                qc.cx(qubit_i,qubit_j)  
                if(dag==True):                                          #                        
                        qc.rz(pi/2,qubit_i)                             #
                        qc.cx(qubit_j,qubit_i)                          #                               n
                        qc.ry(theta/2,qubit_i)                          #  This is the circuit for the B :       
                        qc.cx(qubit_j,qubit_i)                          #                               k
                        qc.ry(-theta/2,qubit_i)                         #
                        qc.rz(-pi/2,qubit_i)                            #
                if(dag==False):                                         #|i>--[X]---Rz(pi/2)Ry(theta/2)-[X]-Ry(-theta/2)--[X]-Rz(-pi/2)---[X]--
                        qc.rz(pi/2,qubit_i)                             #          |                      |                 |            |
                        qc.ry(theta/2,qubit_i)                          #          |                      |                 |            |       
                        qc.cx(qubit_j,qubit_i)                          #|j>------[X]-------------------------------------------------[X]----- 
                        qc.ry(-theta/2,qubit_i)                         # 
                        qc.cx(qubit_j,qubit_i)                          #
                        qc.rz(-pi/2,qubit_i)
                qc.cx(qubit_i,qubit_j)  
                qc.x(qubit_i)
                return qc      
        
        def EigenE(self,e_args):
                n = self.lattice.n
                B = self.B
                half_n = n/2.0
                k = sum(e_args)-half_n
                energy = sqrt((B - cos(2*pi*k/n))**2 + (sin(2*pi*k/n ))**2 )
                return energy
        def evolve_diagonal(self, circuit, time):
                n = self.lattice.n
                B = self.B
                labels= [2,0,1,0]
                circuit.rz(-4*time*B,3)
                for i in range(n-1):
                        i_mode = np.sign(B - cos(2*pi*labels[i]/n))*sqrt((B - cos(2*pi*labels[i]/n))**2 + (sin(2*pi*labels[i]/n ))**2 )
                        circuit.rz(2*time*i_mode,i)
                        circuit.rz(2*time*i_mode,3)
                unitary=print_unitary(circuit,False)
                #print("exp i t U^dag H U")
                #print(unitary)
                #print(circuit.draw())
                return circuit

                        
        def f(self,p,l,k):
                n = self.lattice.n
                e_args = [0]*n
                # l is the length of the original set of indices "s", k is the length of the partion "p" of s           
                for i in p: e_args[i] = 1
                E_k=self.EigenE(e_args) 
                #print("sign = %d l= %d, k = %d"% ((-1)**(l-k),l,k))
                #print(e_args)         
                return (-1)**(l-k)*E_k
        
#==============================================================================================
        #define the coefficient x(i1...in) in the Z-circuit for t-evolution in terms of the eigen-energies
        #the formula for this xs' is given: x(L) = \sum_{K in Comb(L)} (-1)^(l+k) E(K)
        #where l = #L, k = #K, and Comb(X) are all nonempty subsets of X.
#==============================================================================================

        def x(self,s):
                l = len(s)
                n = self.lattice.n
                for k in range(1,l+1):
                        comb = list(itertools.combinations(s, k))       
                        x_s = list(map(lambda elem: self.f([*elem],l,k), comb))
                return np.sum(x_s) 
       

        def Udis(self,qc,dag=False):
                n = self.lattice.n
                B = self.B
                labels= [0,2,1,-1]
                if(dag==True):
                              qc = self.Bogoliubov(qc,2,3,-1,dag)               
                              qc = self.Bogoliubov(qc,0,1,0,dag)               
                              qc = self.Ph_Delay_Beam_Spl(qc,2,3,1, dag)         
                              qc = self.Ph_Delay_Beam_Spl(qc,0,1,0, dag)        
                              qc = self.fSWAP(qc,1,2)                                   
                              qc = self.Ph_Delay_Beam_Spl(qc,2,3,0, dag)        
                              qc = self.Ph_Delay_Beam_Spl(qc,0,1,0, dag)        
                              qc = self.fSWAP(qc,1,2)
                if(dag==False):
                              qc = self.fSWAP(qc,1,2)                   
                              qc = self.Ph_Delay_Beam_Spl(qc,0,1,0, dag)         
                              qc = self.Ph_Delay_Beam_Spl(qc,2,3,0, dag)        
                              qc = self.fSWAP(qc,1,2)                           
                              qc = self.Ph_Delay_Beam_Spl(qc,0,1,0, dag)        
                              qc = self.Ph_Delay_Beam_Spl(qc,2,3,1, dag)        
                              qc = self.Bogoliubov(qc,0,1,0,dag)               
                              qc = self.Bogoliubov(qc,2,3,-1,dag)
                return qc
                
        def ZGateCoefficients(self):
                n = self.lattice.n
                nlist = list(range(n))   
                K_list = []
                for i in range(1,n+1): 
                        comb = list(itertools.combinations(nlist,i))
                        x_comb = list(map(lambda elem: generate_number_op(self.x([*elem]),[*elem],n),comb))                        
                        #x_comb = list(map(lambda elem: (self.x([*elem]),generate_pauli_Z([*elem],n)),comb))                        
                        K_list.append(x_comb)
                        flattened_K_list = list(itertools.chain(*K_list))
                        lattened_K_list = list(itertools.chain(*flattened_K_list))
                flattened_list = [y for x in lattened_K_list for y in x]
                return flattened_list

        def evolve_exact(self,circuit,t,):#coeff_list,barr=False):
                #K_oper = WeightedPauliOperator(coeff_list)
                #circuit = K_oper.evolve(circuit,evo_time=t,num_time_slices=1)
                circuit = self.evolve_diagonal(circuit,t)
                return circuit


# ================================================================ #


class Quench:
      def __init__(self,idx=None,H=None,t0=None,t=None,tstep=None,x0=None):
          self.idx              = idx
          self.H                = H
          self.t0               = t0
          self.t                = t
          self.tstep            = tstep
          self.x0               = x0
          self.backend          = None
          self.quantum_instance = None
      # =====
      def generate_and_run(self,run,print_info=False,shots=1024,layout=[0,1,2,3],optimization=3,calibration=False,output=None):
        if(run=='statevector'):
           # simulazione con algebra lineare e senza barre d'errore
           self.backend = Aer.get_backend('statevector_simulator')
           self.quantum_instance = QuantumInstance(backend=self.backend)
        elif(run=='qasm'):
           # Monte Carlo sulla distribuzione di probabilita' esatta)
           self.backend = Aer.get_backend('qasm_simulator')
           self.quantum_instance = QuantumInstance(backend=self.backend)#,shots=shots)
        else:
           # hardware
           IBMQ.load_account()
           provider = IBMQ.get_provider(hub='ibm-q',group='open',project='main')
           self.backend = provider.get_backend(run)
           if(calibration):
                   self.quantum_instance = QuantumInstance(backend=self.backend,shots=shots,initial_layout=layout,
                                      optimization_level=optimization,measurement_error_mitigation_cls=CompleteMeasFitter)
           else:
               self.quantum_instance = QuantumInstance(backend=self.backend,shots=shots,initial_layout=layout,optimization_level=optimization)
          # =====

        #self.backend          = self.backend
        #self.quantum_instance = self.quantum_instance
        #output.write('field = %f, initial state %s \n'%(self.H.B,self.x0))
        Iave=0

        #----------------------------------------------------------------------------------------------------
	#                                                                 +
	#State |\psi_0> is real time evolved, first  entangling through U      and then  disentangling with U
	#                                                                 dis                          
        #----------------------------------------------------------------------------------------------------

        for time in np.arange(self.t0,self.t,self.tstep):
                qc = generate_initial_circuit(self.x0)
                qc = self.H.Udis(qc,dag = False)
                qc = self.H.evolve_diagonal(qc,time)
                qc = self.H.Udis(qc,dag = True)
                tomography_dict = tomography(output,qc,self.quantum_instance)
                result = process_tomography_dict(output,tomography_dict,self.idx)
                Iave,Istd,S_Aave,S_Astd,S_Bave,S_Bstd,S_ABave,S_ABstd, rho_AB, rho_A, rho_B=result
                eAB,_=LA.eigh(result[8])
                eA,_=LA.eigh(result[9])
                eB,_=LA.eigh(result[10])
                Iave=result[0]
                Istd=result[1]
                output.write("time %f nshots %d I(%s)= %f +/- %f \n" % (time,shots,self.idx,Iave,Istd))
                if(print_info): 
                 ("time %f AB %s \n" % (time,self.idx))
                 output.write('Iave=%f,Istd=%f,S_Aave=%f,S_Astd=%f,S_Bave=%f,S_Bstd=%f,S_ABave=%f,S_ABstd=%f\n, rho_AB=%s,\n rho_A=%s, \n rho_B=%s \n' % result)
                 output.write("rhoAB eigenvalues %s"% eAB)
                 output.write("rhoA eigenvalues %s"% eA)
                 output.write("rhoB eigenvalues %s"% eB)
                output.flush()


import sys
sys.path.append('./')
from sub import *
x0  = [0,0]
c   = 3.0
run_extr= 'statevector'
run_an= 'statevector'
#run = 'qasm'
#quantum_instance_extr = generate_quantum_instance(run_extr,shots=8192,optimization=3,calibration=False)
#quantum_instance_an = generate_quantum_instance(run_an,shots=8192,optimization=3,calibration=False)

logfile_extr = open('sv_out_extr','w')
logfile_extr.write("field= %f initial state = %s hardware = %s"%(c,x0,run_extr))
logfile_an = open('sv_out_an','w')
logfile_an.write("field= %f initial state = %s hardware = %s"%(c,x0,run_an))
for t in np.arange(0,2,0.02): #usa t=6 poi
    #for nshot in [1000,2000,4000,8192]:

        #run = 'qasm'
        quantum_instance_extr = generate_quantum_instance(run_extr,shots=8192,optimization=3,calibration=False)
        quantum_instance_an   = generate_quantum_instance(run_an,  shots=8192,optimization=3,calibration=False)


    # run = 'ibmq_rome' 'ibmq_essex' 'ibmq_ourense' 'ibmq_london' 'ibmq_vigo' 'ibmq_ibmqx2' 'ibmq_burlington' 'ibmq_melbourne'
    # quantum_instance = generate_quantum_instance(run,shots=8192,layout=[*,*],optimization=3,calibration=False)

    # run = 'ibmq_rome' 'ibmq_essex' 'ibmq_ourense' 'ibmq_london' 'ibmq_vigo' 'ibmq_ibmqx2' 'ibmq_burlington' 'ibmq_melbourne'
    # quantum_instance = generate_quantum_instance(run,shots=8192,layout=[*,*],optimization=3,calibration=True)

    # ===========
        #tt=int(t*100)
        #print("time ",tt)
        #print("shot ",nshot)
        qc0 = generate_initial_circuit(x0)
        #qc  = evolve_in_real_time(logfile_extr,qc0,t,c)
        qc1 = exp_itHdis(logfile_an,qc0,t,c)
        #tomography_dict_extr = tomography(logfile_extr,qc, quantum_instance_extr)
        tomography_dict_an   = tomography(logfile_an,qc1,quantum_instance_an)
        #Iave_extr, Istd_extr = process_tomography_dict(logfile_extr,tomography_dict_extr)
        Iave_an, Istd_an     = process_tomography_dict(logfile_an,tomography_dict_an)
        #print('extr I %f %f %f \n' % (t,Iave_extr, Istd_extr) )
        print('an   I %f %f %f \n' % (t,Iave_an, Istd_an) )


        #print("time %f nshot %d I(%s) = %f +/- %f\n"%(t,nshot,x0,Iave,Istd))

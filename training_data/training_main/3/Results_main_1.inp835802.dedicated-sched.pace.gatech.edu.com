from driverConstants import *
from driverStandardMPI import StandardMPIAnalysis
import driverUtils, sys
options = {
    'ams':OFF,
    'analysisType':STANDARD,
    'applicationName':'analysis',
    'aqua':OFF,
    'beamSectGen':OFF,
    'biorid':OFF,
    'cel':OFF,
    'compile_fortran':'ifort -c -fPIC -auto -extend_source -w90 -w95 -WB -I%I -fpp',
    'complexFrequency':OFF,
    'contact':OFF,
    'cosimulation':OFF,
    'coupledProcedure':OFF,
    'cpus':6,
    'cyclicSymmetryModel':OFF,
    'directCyclic':OFF,
    'direct_solver':DMP,
    'double_precision':EXPLICIT,
    'dsa':OFF,
    'dynamic':OFF,
    'filPrt':[],
    'fils':[],
    'finitesliding':OFF,
    'foundation':OFF,
    'geostatic':OFF,
    'heatTransfer':OFF,
    'importer':OFF,
    'importerParts':OFF,
    'includes':('loads.inp',),
    'initialConditionsFile':OFF,
    'input':'main_1',
    'interactive':None,
    'job':'Results_main_1.inp835802.dedicated-sched.pace.gatech.edu',
    'lanczos':OFF,
    'libs':[],
    'massDiffusion':OFF,
    'mp_file_system':(DETECT, DETECT),
    'mp_head_node':('iw-p31-26.pace.gatech.edu', 'iw-p31-26', '172.26.75.50', 'iw-p31-26-ib.pace.gatech.edu', 'iw-p31-26-ib', '172.26.68.111'),
    'mp_host_list':(('iw-p31-26.pace.gatech.edu', 6),),
    'mp_mode':MPI,
    'mp_mode_requested':MPI,
    'mp_mpi_validate':OFF,
    'mp_mpirun_path':'/nv/coraid-pace/usr-local-iw-rhel6/packages/abaqus/6.9/6.9-1/External/mpi/hpmpi-2.2.5.1/bin/mpirun',
    'mp_queue':'PBS',
    'mp_rsh_command':'rsh -n -l pkern3 %H %C',
    'multiphysics':OFF,
    'noDmpDirect':[],
    'noMultiHost':[],
    'no_domain_check':1,
    'outputKeywords':ON,
    'parameterized':OFF,
    'partsAndAssemblies':OFF,
    'parval':OFF,
    'postOutput':OFF,
    'publicSim':OFF,
    'restart':OFF,
    'restartWrite':OFF,
    'rezone':OFF,
    'runCalculator':OFF,
    'scratch':'/nv/hp16/pkern3/scratch/training_data/Results_835802.dedicated-sched.pace.gatech.edu',
    'soils':OFF,
    'solverTypes':['DIRECT', 'DIRECT'],
    'standard_parallel':ALL,
    'staticNonlinear':ON,
    'steadyStateTransport':OFF,
    'step':ON,
    'subGen':OFF,
    'subGenLibs':[],
    'subGenTypes':[],
    'submodel':OFF,
    'substrLibDefs':OFF,
    'substructure':OFF,
    'symmetricModelGeneration':OFF,
    'tmpdir':'/nv/hp16/pkern3/scratch/training_data/Results_835802.dedicated-sched.pace.gatech.edu',
    'tracer':OFF,
    'user':'Al_v144.f',
    'visco':OFF,
}
analysis = StandardMPIAnalysis(options)
status = analysis.run()
sys.exit(status)

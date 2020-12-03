
# Minimal example on how to use the model libraries used in:
#   Lindroos & Hellgren Kotaleski 2020


from   neuron           import h
import MSN_builder          as build
import pickle
import matplotlib.pyplot    as plt
from common_functions import *
# Load model mechanisms
import neuron               as nrn
import numpy as np
nrn.load_mechanisms('mechanisms/single')

h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')

def getParentDendrites(sec, parentsecList):
    if  h.SectionRef(sec=sec).parent and h.SectionRef(sec=sec).parent.name().find('soma') == -1 :
        parentSec = h.SectionRef(sec=sec).parent
        print (" parentSec " + str(parentSec.name()))
        parentsecList.append([parentSec, parentSec.L])

        return getParentDendrites(parentSec, parentsecList)

    else:

        return parentsecList
# specs
specs = {'dspn': {
                    'N': 71,
                    'lib': 'Libraries/D1_71bestFit_updRheob.pkl',
                    'par': 'params_dMSN.json',
                    'morph': 'Morphologies/WT-dMSN_P270-20_1.02_SGA1-m24.swc'},
         'ispn': {
                    'N': 34,
                    'lib': 'Libraries/D2_34bestFit_updRheob.pkl',
                    'par': 'params_iMSN.json',
                    'morph': 'Morphologies/WT-iMSN_P270-09_1.01_SGA2-m1.swc'}
        }

# chose cell type ('ispn' or 'dspn') and model id(s) to simulate
#---------------------------------------------------------------
cell_type         = 'dspn'    # 'dspn'/'ispn'
model_iterator    = range(5)  # range(specs[cell_type]['N']) gives all models

# open library (channel distributions etc)
with open(specs[cell_type]['lib'], 'rb') as f:
    model_sets = pickle.load(f, encoding="latin1")
    print (" model_sets " + str(len(model_sets)))

# simulate model(s)
OUT = {}
for cell_index in model_iterator:

    # initiate cell
    cell = build.MSN(  params=specs[cell_type]['par'],
                       morphology=specs[cell_type]['morph'],
                       variables=model_sets[cell_index]['variables']   )


    # set current injection
    print (model_sets[cell_index])
    rheobase        =   model_sets[cell_index]['rheobase']
    Istim           =   h.IClamp(0.5, sec=cell.soma)
    Istim.delay     =   100
    Istim.dur       =   1000
    Istim.amp       =   (rheobase) *1e-3

    # record vectors
    tm  = h.Vector()
    tm.record(h._ref_t)

    vml = [h.Vector(1e3) for x in cell.dendlist]
    totalCount = 0
    syn_locs = np.array(range(1,9,4))/10 # array([0.1, 0.3, 0.5, 0.7])
    ns      = {}
    nc      = {}
    Syn     = {}
    gbase = 0.3e-3
    delay = 0
    fglut=12.0
    for v,c in zip(vml,cell.dendlist):

        parentDendrites = getParentDendrites(c, [])
        print(" for section " + c.name() + " parents are " + str(parentDendrites) + " dist from soma = " + str(sum([x[1] for x in parentDendrites])))
        # tm[i].record(h._ref_t)
        if totalCount in [21,33]:
            for loc in syn_locs:
            # create a glut synapse (glutamate)
                random_synapse(ns, nc, Syn, c, loc,           \
                                        NS_interval=1000.0/fglut,    \
                                        NC_conductance=gbase,       \
                                        NS_start=delay,             \
                                        seed=None )

        print(c.name())
        v.record(c(0.5)._ref_v)
        totalCount += 1
    # run simulation
    h.finitialize(-80)

    # run simulation
    while h.t < 1000:
        h.fadvance()

    OUT[cell_index] = {}

    plt.clf()
    for v,c in zip(vml,cell.dendlist):
        plt.clf()
        plt.plot(tm.to_python(), v.to_python())
        plt.savefig("output/cell_"+ str(cell_index) + "_dendrite_" + str(c.name()) )

plt.legend()

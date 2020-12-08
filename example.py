
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
import random
nrn.load_mechanisms('mechanisms/single')
h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')

def getParentDendrites(sec, parentsecList):
    if  h.SectionRef(sec=sec).parent and h.SectionRef(sec=sec).parent.name().find('soma') == -1 :
        parentSec = h.SectionRef(sec=sec).parent
        #print (" parentSec " + str(parentSec.name()))
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
model_iterator    = [0]  # range(specs[cell_type]['N']) gives all models

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

    # current needed to make the cell spike (if given to the soma)
    rheobase        =   model_sets[cell_index]['rheobase']
    # Istim           =   h.IClamp(0.5, sec=cell.soma)
    # Istim.delay     =   100
    # Istim.dur       =   1000
    # Istim.amp       =   (rheobase) *1e-3
    # # record vectors
    tm  = h.Vector()
    tm.record(h._ref_t)

    #vml = [h.Vector(1e3) for x in cell.dendlist]
    vm = h.Vector()
    vm.record(cell.soma(0.5)._ref_v)

    # to get count of dendrites
    count_dend = 0
    for x in cell.dendlist:
        count_dend += 1

    rand_dend_nums = random.sample(range(1, len(range(count_dend))), 20)

    print (rand_dend_nums)

    clist = []

    count = 0

    for c in cell.dendlist: # TODO test other dendritic locations!
        print (" count = " + str(count))

        count +=1

        if count not in rand_dend_nums:

            continue

        parentDendrites = getParentDendrites(c, [])
        distFromSoma = sum([x[1] for x in parentDendrites])
        print(parentDendrites)

        count +=1

        # randomly add stim on a distal dendrite
        if len(parentDendrites) > 3 and distFromSoma > 100:

            # create a single glut synapse (glutamate)
            for loc in (np.arange(1,9,.2) )/10:
                syn = h.glutamate( loc, sec=c)
                # default ampa/nmda ratio = 1
                # TODO test other ratios
                syn.ratio = 1
                # TODO vary the time constants
                # tau1_ampa, tau2_ampa, tau1_nmda, tau2_nmda

            # create NetStim object
            # TODO test larger and smaller gbase
            gbase          = 3.0e-3        # 0.3e-3
            ns             = h.NetStim()
            ns.start       = 100
            # TODO test other intervals!
            ns.interval    = 2  # mean interval between two spikes in ms
            ns.noise       = 0  # add noise to ISI? (0-1)
            # TODO change number of epsp
            ns.number      = 10 # numer of sequential psps

            # create NetCon object
            nc             = h.NetCon(ns,syn)
            nc.delay       = 0
            nc.weight[0]   = gbase

            # record local potential
            vml = h.Vector() # record local potential
            vml.record(c(0.5)._ref_v)

            # record synpatic currents # TODO record conductance
            ampa_current = h.Vector()
            ampa_current.record(syn._ref_i_ampa)
            nmda_current = h.Vector()
            nmda_current.record(syn._ref_i_nmda)

            clist.append(c)
            # break   # set in first only...
            # if count > 4:
            #     break
    # run simulation
    h.finitialize(-85)

    # run simulation
    while h.t < 300:
        h.fadvance()

    print (clist)

    for c in clist:
        plt.clf()
        print (" plotting for c = " + c.name())
        fig,ax = plt.subplots(2,1, figsize=(6,10), sharex='all')
        ax[0].plot(tm.to_python(), vm.to_python(), label='soma')
        ax[0].plot(tm.to_python(), vml.to_python(), label=c.name())
        ax[1].plot(tm, ampa_current, label='ampa')
        ax[1].plot(tm, nmda_current, label='nmda')

        ax[0].set_title('Voltage in the soma and a stimulated dendrite')
        ax[0].set_ylabel('Voltage (mV)')
        ax[0].legend()
        ax[1].set_title('synaptic current')
        ax[1].set_ylabel('Current (nA)')
        ax[1].set_xlabel('Time (ms)')
        ax[1].legend()

        plt.savefig("output/final_out_" + c.name() + ".png")

    #plt.savefig("output/cell_"+ str(cell_index) + "_dendrite_" + str(c.name()) )

plt.legend()

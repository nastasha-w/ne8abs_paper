
import numpy as np
import os

from fire_an.ionrad.ion_utils import Linetable_PS20
import fire_an.mainfunc.makehist as mh
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c
    
def run_hist_ptmasses_all2(opt):
    if opt >= 0 and opt < 450:
        # z = 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
        # m13-sr
        # 450 indices
        ind = opt - 0
        simnames = sl.m13_sr_all2 # len 15
        snaps = sl.snaps_sr # len 6
        pts = [0, 1, 2, 4, 5] # len 5
    elif opt >= 450 and opt < 510:
        # z = 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
        # m13-hr
        # 60 indices
        ind = opt - 450
        simnames = sl.m13_hr_all2 # len 2
        snaps = sl.snaps_hr # len 6
        pts = [0, 1, 2, 4, 5] # len 5
    elif opt >= 510 and opt < 630:
        # z = 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
        # m12-sr
        # 120 indices
        ind = opt - 510
        simnames = sl.m12_sr_all2 # len 4
        snaps = sl.snaps_sr # len 6
        pts = [0, 1, 2, 4, 5] # len 5
    elif opt >= 630 and opt < 1170:
        # z = 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
        # m12-hr
        # 540 indices
        ind = opt - 630
        simnames = sl.m12_hr_all2 # len 18
        snaps = sl.snaps_hr # len 6
        pts = [0, 1, 2, 4, 5] # len 5
    
    _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
    outdir = '/scratch1/08466/tg877653/output/hists/ptmasses_all2/'
    outname = 'hist_r3D_by_mass_pt{pt}_{simname}_snap{snap}_bins1_v1.hdf5'
    simi = ind // (len(snaps) * len(pts))
    snpi = (ind % (len(snaps) * len(pts))) // (len(pts))
    pti = (ind % len(pts))
    simname = simnames[simi]
    snapnum = snaps[snpi]
    pt = pts[pti]

    if pt == 5 and 'sdp1e10' in simname:
        msg = (f'Skipping particle type {pt} for noBH'
                f' simulation {simname}, snap {snapnum}')
        print(msg)
        return None
    
    runit = 'Rvir'
    rbins = np.arange(0.15, 5., 0.01)
    rbins = np.append(np.arange(0., 0.15, 0.005), rbins)
    weighttype = 'Mass' 
    weighttype_args = {}
    axtypes = []
    axtypes_args = []
    axbins = []

    # directory is halo name + resolution 
    dp2 = '_'.join(simname.split('_')[:2])
    if dp2.startswith('m13h02_'):
        dp2 = dp2.replace('m13h02', 'm13h002')
    dirpath = '/'.join([_dirpath, dp2, simname])

    outfilen = outdir + outname.format(pt=pt, simname=simname, 
                                       snap=snapnum)
    mh.histogram_radprof(dirpath, snapnum,
                         weighttype, weighttype_args, axtypes, axtypes_args,
                         particle_type=pt, 
                         center='shrinksph', rbins=rbins, runit=runit,
                         logweights=True, logaxes=True, axbins=axbins,
                         outfilen=outfilen, overwrite=False)
    return None

def run_hist_o6ne8mg10_tosample2(opt):
    outdir = '/scratch1/08466/tg877653/output/hists/all2_model3/'
    outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
    particle_type = 0
    _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
    if opt >= 0 and opt < 1824:
        # 1824 inds, 864 with ion weights
        # m12-hr (m12f already run; complete to all2 set)
        iw_base = 0
        niw_base = 864
        snaps = sl.snaps_hr # z=0.5 - 1.0, len 6
        # len 16
        simnames = [('m12b_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp2e-4_gacc31_fa0.5'),
                    ('m12c_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e-4_gacc31_fa0.5'),
                    ('m12i_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp2e-4_gacc31_fa0.5'),
                    ('m12m_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp2e-4_gacc31_fa0.5'),
                    ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp2e-4_gacc31_fa0.5'),
                    ('m12r_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e-4_gacc31_fa0.5'),
                    ('m12w_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e-4_gacc31_fa0.5'),
                    ('m12z_m4e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e-4_gacc31_fa0.5'),
                    ('m12b_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m12c_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m12i_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m12m_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m12r_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m12w_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m12z_m4e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ]
    elif opt >= 1824 and opt < 2166:
        # 342 inds, 162 with ion weights
        # m12-sr (m12f already run; complete to all2 set)
        iw_base = 1824
        niw_base = 1986
        snaps = sl.snaps_sr # z=0.5 - 1.0, len 6
        # len 3
        simnames = [('m12i_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                    ('m12m_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                    ('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                    ]
    # all m13-hr (h113, h206) already run
    elif opt >= 2166 and opt < 3420:
        # 1254 inds, 594 with ion weights
        # m13-sr (m13h113, m13h206 already run; complete to all2 set)
        iw_base = 2166
        niw_base = 2760
        snaps = sl.snaps_sr # z=0.5 - 1.0, len 6
        # len 11
        simnames = [('m13h002_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m13h007_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m13h029_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m13h223_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m13h236_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m13h002_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                    ('m13h007_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                    ('m13h009_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'), 
                    ('m13h029_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                    ('m13h037_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'), 
                    ('m13h236_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                    ]
    # total: len 19
    if (opt >= iw_base and opt < niw_base): # len 9
        ind = opt - iw_base
        wts = ['O6', 'Ne8', 'Mg10'] # len 3
        axtypes_opts = [['sim-direct']] * 3 # len 3
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/{elt}'}],
                            ]
        axqts = ['Temperature', 'Density', '{elt}']
        axbins = [0.05, 0.05, 0.1]
    else: # len 10
        ind = opt - niw_base
        wts = ['Mass', 'Volume'] # len 2
        axtypes_opts = [['sim-direct']] * 5 # len 5
        axtypes_args_opts = [[{'field': 'Temperature'}],
                                [{'field': 'Density'}],
                                [{'field': 'ElementAbundance/Oxygen'}],
                                [{'field': 'ElementAbundance/Neon'}],
                                [{'field': 'ElementAbundance/Magnesium'}],
                            ]
        axqts = ['Temperature', 'Density', 'Oxygen', 'Neon', 'Magnesium']
        axbins = [0.05, 0.05] + [0.1] * 3

    simi = ind // (len(snaps) * len(wts) * len(axqts))
    snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
            // (len(wts) * len(axqts))
    wti = (ind % (len(wts) * len(axqts))) // len(axqts)
    axi = ind % len(axqts)
    simname = simnames[simi]
    snapnum = snaps[snpi]
    wt = wts[wti]
    axtypes = axtypes_opts[axi]
    axtypes_args = axtypes_args_opts[axi]
    axqt = axqts[axi]
    
    runit = 'pkpc'
    rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
            np.arange(40., 1001., 20.)
    rbins = np.append(np.arange(0., 40., 5.), rbins)

    # directory is halo name + resolution 
    dp2 = '_'.join(simname.split('_')[:2])
    if dp2.startswith('m13h02_'):
        dp2 = dp2.replace('m13h02', 'm13h002')
    dirpath = '/'.join([_dirpath, dp2, simname])

    if wt in ['Mass', 'Volume']:
        weighttype = wt
        weighttype_args = dict()
    else:
        weighttype = 'ion'
        weighttype_args = {'ps20depletion': False, 'ion': wt,
                           'density': False}
        if wt == 'H1':
            weighttype_args.update({'ionfrac-method': 'sim'})
        else:
            dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                        vol=True, lintable=True)
            parentelt = dummytab.element
            axtypes_args = \
                [{key: (dct[key]).format(elt=parentelt) for key in dct}
                  for dct in axtypes_args]
            axqt = axqt.format(elt=parentelt)
    outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                       snap=snapnum)

    mh.histogram_radprof(dirpath, snapnum,
                         weighttype, weighttype_args, axtypes, axtypes_args,
                         particle_type=particle_type, 
                         center='shrinksph', rbins=rbins, runit=runit,
                         logweights=True, logaxes=True, axbins=axbins,
                         outfilen=outfilen, overwrite=False)


def run_phasediagrams_radius(opt):
    # nH-T weighted by M, V, Ne, Ne8
    # sample all2
    # 4 weights -> 4 runs per sim/snap
    outdir = '/scratch1/08466/tg877653/output/hists/phasediagrams_all2/'
    
    if opt >= 0  and opt < 240:
        # outdir stampede2
        outdir = '/scratch/08466/tg877653/output/hists/phasediagrams_all2/'
        # 240 indices
        ind = opt - 0
        simnames = sl.m12_f2md # len 8, + 2 for crheatfix
        snaps = sl.snaps_f2md # len 6
    elif opt >= 240 and opt < 336:
        # 96 indices
        ind = opt - 240
        simnames = sl.m12_sr_all2 # len 4
        snaps = sl.snaps_sr # len 6
    elif opt >= 336 and opt < 768:
        # 432 indices
        ind = opt - 336
        simnames = sl.m12_hr_all2 # len 18
        snaps = sl.snaps_hr # len 6
    elif opt >= 768 and opt < 1128:
        # 360 indices
        ind = opt - 768
        simnames = sl.m13_sr_all2 # len 15
        snaps = sl.snaps_sr # len 6
    elif opt >= 1128 and opt < 1176:
        # 48 indices
        ind = opt - 1128
        simnames = sl.m13_hr_all2 # len 2
        snaps = sl.snaps_hr # len 6
    
    wts = ['Mass', 'Volume', 'Metal'] + ['ion'] 
    wtargs = [{}, {}, 
              {'element': 'Neon'},
              {'ion': 'Ne8', 'ps20depletion': False},
              ]
    runit = 'Rvir'
    rbins = np.linspace(0.0, 1.3, 27)
    at = ['sim-direct', 'Metal']
    atarg = [{'field': 'Temperature'},
             {'element': 'Hydrogen', 'density': True}]
    axbin = [0.05, 0.05]
    logaxes = [True, True]

    #_dirpath = '/scratch/projects/xsede/GalaxiesOnFIRE/metal_diffusion/'
    simi = ind // (len(snaps) * len(wts))
    snpi = (ind % (len(snaps) * len(wts))) \
           // (len(wts))
    wti = ind % (len(wts))

    simname = simnames[simi]
    snapnum = snaps[snpi]
    wt = wts[wti]
    wtarg = wtargs[wti]

    #dirpath = '/'.join([_dirpath, simname])
    dirpath = sl.dirpath_from_simname(simname)

    atstr = 'rcen_temperature_density'
    wtstr = 'gasmass' if wt == 'Mass' else\
            'gasvol' if wt == 'Volume' else\
            wtarg['ion'] if wt == 'ion' else \
            wtarg['element']
    outfilen = outdir +\
               (f'hist_{atstr}_by_{wtstr}_{simname}_snap{snapnum}'
                '_bins1_v1.hdf5')
    if os.path.isfile(outfilen):
        print(outfilen, ' already exists; skipping')
        return None
    mh.histogram_radprof(dirpath, snapnum,
                         wt, wtarg, at, atarg,
                         particle_type=0, 
                         center='shrinksph', rbins=rbins, runit=runit,
                         logweights=True, logaxes=logaxes, axbins=axbin,
                         outfilen=outfilen, overwrite=False)

def run_hist_Zprof(opt):
    # for total halo metallicities
    # 1 run per sim/snap  
    outdir = '/scratch1/08466/tg877653/output/hists/r_wtd/'
    if opt >= 0 and opt < 60:
        # 60 indices
        outdir = '/projects/b1026/nastasha/hists/r_wtd/'
        ind = opt - 0
        simnames = sl.m12_f2md # len 8, + 2 for crheatfix
        snaps = sl.snaps_f2md # len 6
    elif opt >= 60 and opt < 84:
        # 24 indices
        ind = opt - 60
        simnames = sl.m12_sr_all2 # len 4
        snaps = sl.snaps_sr # len 6
    elif opt >= 84 and opt < 192:
        # 108 indices
        ind = opt - 84
        simnames = sl.m12_hr_all2 # len 18
        snaps = sl.snaps_hr # len 6
    elif opt >= 192 and opt < 282:
        # 90 indices
        ind = opt - 192
        simnames = sl.m13_sr_all2 # len 15
        snaps = sl.snaps_sr # len 6
    elif opt >= 282 and opt < 294:
        # 12 indices
        ind = opt - 282
        simnames = sl.m13_hr_all2 # len 2
        snaps = sl.snaps_hr # len 6
    elif opt >= 294 and opt < 312:
        outdir = '/projects/b1026/nastasha/hists/r_wtd/'
        # 18 indices
        ind = opt - 294
        simnames = sl.m12plus_f3nobh
        snaps = sl.snaps_hr
    elif opt >= 312 and opt < 366:
        outdir = '/projects/b1026/nastasha/hists/r_wtd/'
        # 54 indices
        ind = opt - 312
        simnames = sl.m12plus_f3nobh_lores
        snaps = sl.snaps_hr
    elif opt >= 366 and opt < 390:
        # 24 indices
        ind = opt - 366
        simnames = sl.m12_fire3x_tests
        snaps = sl.snaps_sr

    #_dirpath = '/scratch/projects/xsede/GalaxiesOnFIRE/metal_diffusion/'
    simi = ind // (len(snaps))
    snpi = ind % len(snaps)
    simname = simnames[simi]
    snapnum = snaps[snpi]

    at = ['sim-direct', 'sim-direct']
    atarg = [{'field': 'Metallicity'}, {'field': 'Temperature'}]
    axbin = [0.05, 0.1]
    logaxes = [True, True]
    wt = 'Mass'
    wtarg = {}

    runit = 'Rvir'
    rbins = np.linspace(0.0, 1.3, 27)

    #dirpath = '/'.join([_dirpath, simname])
    dirpath = sl.dirpath_from_simname(simname)

    atstr = 'rcen_Metallicity_Temperature'
    wtstr = 'gasmass'
    outfilen = outdir +\
               (f'hist_{atstr}_by_{wtstr}_{simname}_snap{snapnum}'
                '_bins1_v1_hvcen.hdf5')
    if os.path.isfile(outfilen):
        print(outfilen, ' already exists; skipping')
        return None
    mh.histogram_radprof(dirpath, snapnum,
                         wt, wtarg, at, atarg,
                         particle_type=0, 
                         center='shrinksph', rbins=rbins, runit=runit,
                         logweights=True, logaxes=logaxes, axbins=axbin,
                         outfilen=outfilen, overwrite=False)

def run_hist_mstellar_Zstellar(opt):
    # for total halo stellar masses
    # 1 run per sim/snap  
    outdir = '/scratch1/08466/tg877653/output/hists/r_wtd/'
    if opt >= 0 and opt < 60:
        # 60 indices
        outdir = '/projects/b1026/nastasha/hists/r_wtd/'
        ind = opt - 0
        simnames = sl.m12_f2md # len 8, + 2 for crheatfix
        snaps = sl.snaps_f2md # len 6
    elif opt >= 60 and opt < 84:
        # 24 indices
        ind = opt - 60
        simnames = sl.m12_sr_all2 # len 4
        snaps = sl.snaps_sr # len 6
    elif opt >= 84 and opt < 192:
        # 108 indices
        ind = opt - 84
        simnames = sl.m12_hr_all2 # len 18
        snaps = sl.snaps_hr # len 6
    elif opt >= 192 and opt < 282:
        # 90 indices
        ind = opt - 192
        simnames = sl.m13_sr_all2 # len 15
        snaps = sl.snaps_sr # len 6
    elif opt >= 282 and opt < 294:
        # 12 indices
        ind = opt - 282
        simnames = sl.m13_hr_all2 # len 2
        snaps = sl.snaps_hr # len 6
    elif opt >= 294 and opt < 312:
        outdir = '/projects/b1026/nastasha/hists/r_wtd/'
        # 18 indices
        ind = opt - 294
        simnames = sl.m12plus_f3nobh
        snaps = sl.snaps_hr
    elif opt >= 312 and opt < 366:
        outdir = '/projects/b1026/nastasha/hists/r_wtd/'
        # 54 indices
        ind = opt - 312
        simnames = sl.m12plus_f3nobh_lores
        snaps = sl.snaps_hr
    elif opt >= 366 and opt < 390:
        # 24 indices
        ind = opt - 366
        simnames = sl.m12_fire3x_tests
        snaps = sl.snaps_sr

    #_dirpath = '/scratch/projects/xsede/GalaxiesOnFIRE/metal_diffusion/'
    simi = ind // (len(snaps))
    snpi = ind % len(snaps)
    simname = simnames[simi]
    snapnum = snaps[snpi]

    at = ['sim-direct']
    atarg = [{'field': 'Metallicity'}]
    axbin = [0.05]
    logaxes = [True]
    wt = 'Mass'
    wtarg = {}

    runit = 'Rvir'
    rbins = np.linspace(0.0, 1.3, 53)

    #dirpath = '/'.join([_dirpath, simname])
    dirpath = sl.dirpath_from_simname(simname)

    atstr = 'rcen_'
    wtstr = 'stellarmass'
    outfilen = outdir +\
               (f'hist_{atstr}_by_{wtstr}_{simname}_snap{snapnum}'
                '_bins1_v1_hvcen.hdf5')
    if os.path.isfile(outfilen):
        print(outfilen, ' already exists; skipping')
        return None
    mh.histogram_radprof(dirpath, snapnum,
                         wt, wtarg, at, atarg,
                         particle_type=4, 
                         center='shrinksph', rbins=rbins, runit=runit,
                         logweights=True, logaxes=logaxes, axbins=axbin,
                         outfilen=outfilen, overwrite=False)

def run_hist(opt):
    if opt >= 0 and opt < 60:
        ind = opt
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['Mass', 'Volume', 'H1']
        snaps = [50] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h113_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    'm13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                   ]
        axtypes_opts = [['sim-direct']] * 5
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/Oxygen'}],
                             [{'field': 'ElementAbundance/Neon'}],
                             [{'field': 'ElementAbundance/Magnesium'}],
                            ]
        axqts = ['Temperature', 'Density', 'Oxygen', 'Neon', 'Magnesium']
        axbins = [0.05, 0.05] + [0.1] * 3

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 60 and opt < 96:
        ind = opt - 60
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['O6', 'Ne8', 'Mg10']
        snaps = [50] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h113_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    'm13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                   ]
        axtypes_opts = [['sim-direct']] * 3
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/{elt}'}],
                            ]
        axqts = ['Temperature', 'Density', '{elt}']
        axbins = [0.05, 0.05, 0.1]

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)

        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 96 and opt < 126:
        ind = opt - 96
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['Mass', 'Volume', 'H1']
        snaps = [258] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
                    'm12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
                    ]
        axtypes_opts = [['sim-direct']] * 5
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/Oxygen'}],
                             [{'field': 'ElementAbundance/Neon'}],
                             [{'field': 'ElementAbundance/Magnesium'}],
                            ]
        axqts = ['Temperature', 'Density', 'Oxygen', 'Neon', 'Magnesium']
        axbins = [0.05, 0.05] + [0.1] * 3

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 126 and opt < 144:
        ind = opt - 126
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['O6', 'Ne8', 'Mg10']
        snaps = [258] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
                    'm12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
                    ]
        axtypes_opts = [['sim-direct']] * 3
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/{elt}'}],
                            ]
        axqts = ['Temperature', 'Density', '{elt}']
        axbins = [0.05, 0.05, 0.1]

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)

        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 144 and opt < 174:
        ind = opt - 144
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['Mass', 'Volume', 'H1']
        snaps = [258] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h113_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31_fa0.5',
                    'm13h206_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp3e-4_gacc31_fa0.5',
                   ]
        axtypes_opts = [['sim-direct']] * 5
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/Oxygen'}],
                             [{'field': 'ElementAbundance/Neon'}],
                             [{'field': 'ElementAbundance/Magnesium'}],
                            ]
        axqts = ['Temperature', 'Density', 'Oxygen', 'Neon', 'Magnesium']
        axbins = [0.05, 0.05] + [0.1] * 3

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 174 and opt < 192:
        ind = opt - 174
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['O6', 'Ne8', 'Mg10']
        snaps = [258] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h113_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31_fa0.5',
                    'm13h206_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp3e-4_gacc31_fa0.5',
                   ]
        axtypes_opts = [['sim-direct']] * 3
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/{elt}'}],
                            ]
        axqts = ['Temperature', 'Density', '{elt}']
        axbins = [0.05, 0.05, 0.1]

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    
    elif opt >= 192 and opt < 207:
        ind = opt - 192
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['Mass', 'Volume', 'H1']
        snaps = [50] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    ]
        axtypes_opts = [['sim-direct']] * 5
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/Oxygen'}],
                             [{'field': 'ElementAbundance/Neon'}],
                             [{'field': 'ElementAbundance/Magnesium'}],
                            ]
        axqts = ['Temperature', 'Density', 'Oxygen', 'Neon', 'Magnesium']
        axbins = [0.05, 0.05] + [0.1] * 3

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 207 and opt < 216:
        ind = opt - 206
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['O6', 'Ne8', 'Mg10']
        snaps = [50] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    ]
        axtypes_opts = [['sim-direct']] * 3
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/{elt}'}],
                            ]
        axqts = ['Temperature', 'Density', '{elt}']
        axbins = [0.05, 0.05, 0.1]

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    
    # add more redshifts (set 2) to previous sample
    #----------------------------------------------
    elif opt >= 216 and opt < 516:
        ind = opt - 216
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['Mass', 'Volume', 'H1']
        snaps = [45, 46, 47, 48, 49] # z=1.0 - 0.6
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h113_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    'm13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                   ]
        axtypes_opts = [['sim-direct']] * 5
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/Oxygen'}],
                             [{'field': 'ElementAbundance/Neon'}],
                             [{'field': 'ElementAbundance/Magnesium'}],
                            ]
        axqts = ['Temperature', 'Density', 'Oxygen', 'Neon', 'Magnesium']
        axbins = [0.05, 0.05] + [0.1] * 3

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 516 and opt < 696:
        ind = opt - 516
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['O6', 'Ne8', 'Mg10']
        snaps = [45, 46, 47, 48, 49] # z=1.0 - 0.6
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h113_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    'm13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                   ]
        axtypes_opts = [['sim-direct']] * 3
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/{elt}'}],
                            ]
        axqts = ['Temperature', 'Density', '{elt}']
        axbins = [0.05, 0.05, 0.1]

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)

        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 696 and opt < 846:
        ind = opt - 696
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['Mass', 'Volume', 'H1']
        snaps = [186, 197, 210, 224, 240] # z=1.0 - 0.6
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
                    'm12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
                    ]
        axtypes_opts = [['sim-direct']] * 5
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/Oxygen'}],
                             [{'field': 'ElementAbundance/Neon'}],
                             [{'field': 'ElementAbundance/Magnesium'}],
                            ]
        axqts = ['Temperature', 'Density', 'Oxygen', 'Neon', 'Magnesium']
        axbins = [0.05, 0.05] + [0.1] * 3

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 846 and opt < 936:
        ind = opt - 846
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['O6', 'Ne8', 'Mg10']
        snaps = [186, 197, 210, 224, 240] # z=1.0 - 0.6
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
                    'm12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
                    ]
        axtypes_opts = [['sim-direct']] * 3
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/{elt}'}],
                            ]
        axqts = ['Temperature', 'Density', '{elt}']
        axbins = [0.05, 0.05, 0.1]

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)

        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 936 and opt < 1086:
        ind = opt - 936
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['Mass', 'Volume', 'H1']
        snaps = [186, 197, 210, 224, 240] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h113_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31_fa0.5',
                    'm13h206_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp3e-4_gacc31_fa0.5',
                   ]
        axtypes_opts = [['sim-direct']] * 5
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/Oxygen'}],
                             [{'field': 'ElementAbundance/Neon'}],
                             [{'field': 'ElementAbundance/Magnesium'}],
                            ]
        axqts = ['Temperature', 'Density', 'Oxygen', 'Neon', 'Magnesium']
        axbins = [0.05, 0.05] + [0.1] * 3

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 1086 and opt < 1176:
        ind = opt - 1086
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['O6', 'Ne8', 'Mg10']
        snaps = [186, 197, 210, 224, 240] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h113_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31_fa0.5',
                    'm13h206_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp3e-4_gacc31_fa0.5',
                   ]
        axtypes_opts = [['sim-direct']] * 3
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/{elt}'}],
                            ]
        axqts = ['Temperature', 'Density', '{elt}']
        axbins = [0.05, 0.05, 0.1]

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    
    elif opt >= 1176 and opt < 1251:
        ind = opt - 1176
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['Mass', 'Volume', 'H1']
        snaps = [45, 46, 47, 48, 49] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    ]
        axtypes_opts = [['sim-direct']] * 5
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/Oxygen'}],
                             [{'field': 'ElementAbundance/Neon'}],
                             [{'field': 'ElementAbundance/Magnesium'}],
                            ]
        axqts = ['Temperature', 'Density', 'Oxygen', 'Neon', 'Magnesium']
        axbins = [0.05, 0.05] + [0.1] * 3

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 1251 and opt < 1296:
        ind = opt - 1251
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['O6', 'Ne8', 'Mg10']
        snaps = [45, 46, 47, 48, 49] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    ]
        axtypes_opts = [['sim-direct']] * 3
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/{elt}'}],
                            ]
        axqts = ['Temperature', 'Density', '{elt}']
        axbins = [0.05, 0.05, 0.1]

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)

    mh.histogram_radprof(dirpath, snapnum,
                         weighttype, weighttype_args, axtypes, axtypes_args,
                         particle_type=particle_type, 
                         center='shrinksph', rbins=rbins, runit=runit,
                         logweights=True, logaxes=True, axbins=axbins,
                         outfilen=outfilen, overwrite=False)

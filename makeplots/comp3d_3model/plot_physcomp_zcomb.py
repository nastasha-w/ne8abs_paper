'''
Each panel shows redshift ensemble medians and scatter for different
physics models. The rows show the properties density, temperature, 
and metallicity. Different column show different weights.
There are different plots for different weight/property sets and ICs.
'''
import h5py
import numpy as np

import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt

import ne8abs_paper.makeplots.plot_utils as pu
import ne8abs_paper.makeplots.tol_colors as tc
import ne8abs_paper.utils.constants_and_units as c
import ne8abs_paper.utils.opts_locs as ol

def plot3dprof_physcomp_weightvalgrid(filen_template, fillkw_weights, 
                                      fillkw_vals, fillkw_z, fillkw_phys,
                                      figtitle=None, outname=None):
    rholab = ('$\\log_{10}\\,\\rho \\; '
              '[\\mathrm{H}(X=0.752) \\, \\mathrm{cm}^{-3}]$')
    tlab = '$\\log_{10}\\,\\mathrm{T} \\; [\\mathrm{K}]$'
    zolab = ('$\\log_{10}\\,\\mathrm{Z}_\\mathrm{O} \\;'
             ' [\\mathrm{Z}_{\\mathrm{O}, \\odot}]$')
    znelab = ('$\\log_{10}\\,\\mathrm{Z}_\\mathrm{Ne} \\;'
              ' [\\mathrm{Z}_{\\mathrm{Ne}, \\odot}]$')
    zmglab = ('$\\log_{10}\\,\\mathrm{Z}_\\mathrm{Mg} \\;'
              ' [\\mathrm{Z}_{\\mathrm{Mg}, \\odot}]$')
    ylabopts = {('sim-direct', 'Density'): rholab,
                ('sim-direct', 'Temperature'): tlab,
                ('sim-direct', 'ElementAbundance/Oxygen'): zolab,
                ('sim-direct', 'ElementAbundance/Neon'): znelab,
                ('sim-direct', 'ElementAbundance/Magnesium'): zmglab,
               }
    rhonorm = c.atomw_H * c.u / 0.752
    tnorm = 1.
    zonorm = ol.solar_abunds_ea['oxygen']
    znenorm = ol.solar_abunds_ea['neon']
    zmgnorm = ol.solar_abunds_ea['magnesium']
    ynormopts = {('sim-direct', 'Density'): rhonorm,
                 ('sim-direct', 'Temperature'): tnorm,
                 ('sim-direct', 'ElementAbundance/Oxygen'): zonorm,
                 ('sim-direct', 'ElementAbundance/Neon'): znenorm,
                 ('sim-direct', 'ElementAbundance/Magnesium'): zmgnorm,
                }
    wlabopts = {'Mass': 'Mass-wtd',
                'Volume': 'Volume-wtd',
                ('ion', 'O6'): 'O VI-wtd',
                ('ion', 'Ne8'): 'Ne VIII-wtd',
                ('ion', 'Mg10'): 'Mg X-wtd',
                ('ion', 'H1'): 'H I-wtd',
                }
    cset = tc.tol_cset('bright')
    pcolors = {'AGN-CR': cset.green,
               'AGN-noCR': cset.red,
               'noBH': cset.blue,
               }
    hists_raw = []
    rvals = []
    yvals = []
    ykeys = []
    wkeys = []
    pkeys = []
    rvirs = []
    zvals = []
    runits = None
    for vdict in fillkw_vals:
        tdict = vdict.copy()
        _rvals = []
        _yvals = []
        _ykeys = []
        _wkeys = []
        _pkeys = []
        _rvirs = []
        _zvals = []
        _hists_raw = []
        for wdict in fillkw_weights:
            _tdict = tdict.copy()
            _tdict.update(wdict)
            __rvals = []
            __yvals = []
            __ykeys = []
            __wkeys = []
            __pkeys = []
            __rvirs = []
            __zvals = []
            __hists_raw = []
            for pdict, zlist in zip(fillkw_phys, fillkw_z):
                __tdict = _tdict.copy()
                __tdict.update(pdict)
                ___rvals = []
                ___yvals = []
                ___ykeys = []
                ___wkeys = []
                ___pkeys = []
                ___rvirs = []
                ___zvals = []
                ___hists_raw = []
                for zdict in zlist:
                    ___tdict = __tdict.copy()
                    ___tdict.update(zdict)
                    filen = filen_template.format(**___tdict)
                    with h5py.File(filen, 'r') as f:
                        hdatpath = 'Header/inputpars/halodata'
                        rvir_cm = f[hdatpath].attrs['Rvir_cm']
                        rv = f['axis_0/bins'][:]
                        _runits = f['axis_0'].attrs['units'].decode()
                        yv = f['axis_1/bins'][:]
                        ykey = f['axis_1'].attrs['qty_type'].decode()
                        if ykey == 'sim-direct':
                            _path = 'axis_1/qty_type_args'
                            _key = 'field'
                            ykey = (ykey, f[_path].attrs[_key].decode())
                        elif ykey == 'ion':
                            _path = 'axis_1/qty_type_args'
                            _key = 'ion'
                            ykey = (ykey, f[_path].attrs[_key].decode())
                        elif ykey == 'Metal':
                            _path = 'axis_1/qty_type_args'
                            _key = 'element'
                            ykey = (ykey, f[_path].attrs[_key].decode())
                        hist_raw = f['histogram/histogram'][:]
                        wkey = f['histogram'].attrs['weight_type'].decode()
                        if wkey == 'ion':
                            _path = 'histogram/weight_type_args'
                            _key = 'ion'
                            wkey = (wkey, f[_path].attrs[_key].decode())
                        elif wkey == 'Metal':
                            _path = 'histogram/weight_type_args'
                            _key = 'element'
                            wkey = (wkey, f[_path].attrs[_key].decode())
                        pkey = 'noBH' if '_sdp1e10_' in filen else \
                           'AGN-CR' if '_MHDCRspec1_' in filen else \
                           'AGN-noCR' if '_MHD_' in filen else \
                            None
                        redshift = f['Header/cosmopars'].attrs['z']
                        
                        # adjust for plotting purposes
                        if yv[0] == -np.inf:
                            yv[0] = 2. * yv[1] - yv[2]
                        if yv[-1] == np.inf:
                            yv[-1] = 2. * yv[-2] - yv[-3]
                        
                        if runits is not None:
                            if _runits != runits:
                                raise ValueError('mismatch in r3D units')
                        if _runits == 'pkpc':
                            rvir = rvir_cm / (c.cm_per_mpc * 1e-3)
                        elif _runits == 'Rvir':
                            rvir = 1.
                    ___hists_raw.append(hist_raw)
                    ___rvals.append(rv)
                    ___yvals.append(yv)
                    ___ykeys.append(ykey)
                    ___wkeys.append(wkey)
                    ___pkeys.append(pkey)
                    ___rvirs.append(rvir)
                    ___zvals.append(redshift)
                __hists_raw.append(___hists_raw)
                __rvals.append(___rvals)
                __yvals.append(___yvals)
                __ykeys.append(___ykeys)
                __wkeys.append(___wkeys)
                __pkeys.append(___pkeys)
                __rvirs.append(___rvirs)
                __zvals.append(___zvals)
            _hists_raw.append(__hists_raw)
            _rvals.append(__rvals)
            _yvals.append(__yvals)
            _ykeys.append(__ykeys)
            _wkeys.append(__wkeys)
            _pkeys.append(__pkeys)
            _rvirs.append(__rvirs)
            _zvals.append(__zvals)
        hists_raw.append(_hists_raw)
        rvals.append(_rvals)
        yvals.append(_yvals)
        ykeys.append(_ykeys)
        wkeys.append(_wkeys)
        pkeys.append(_pkeys)
        rvirs.append(_rvirs)
        zvals.append(_zvals)
    #print(ykeys)
    #print(wkeys)

    if runits == 'Rvir':
        xlabel = '$\\mathrm{r} \\; [\\mathrm{R}_{\\mathrm{vir}}]$'
        markrvir = False
    else:
        xlabel = '$\\mathrm{r} \\; [\\mathrm{pkpc}]$'
        markrvir = True

    panelsize = 2.
    legheight = 2.
    ncols = len(fillkw_weights)
    nrows = len(fillkw_vals)
    wspace = 0.
    hspace = 0. 
    nph = len(fillkw_phys)
    nz = len(fillkw_z[0])
    width_ratios = [panelsize] * ncols
    height_ratios = [panelsize] * nrows + [legheight]
    width = sum(width_ratios) \
            * (1. + (len(width_ratios) - 1.) / len(width_ratios) * wspace)
    height = sum(height_ratios) \
             * (1. + (len(height_ratios) - 1.) / len(height_ratios) * hspace)

    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(nrows=nrows + 1, ncols=ncols, hspace=hspace, 
                        wspace=wspace, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [[fig.add_subplot(grid[i, j]) \
             for j in range(ncols)] for i in range(nrows)]
    lax = fig.add_subplot(grid[nrows, :])
    fontsize = 12

    if not np.all([np.all([np.all([np.all([key == rowkeys[0][0][0] \
                                           for key in _key_list])\
                                   for _key_list in key_list])\
                           for key_list in rowkeys])\
                   for rowkeys in ykeys]):
        print('ykeys:\n', ykeys)
        raise RuntimeError('y-axis quantities differ in same row')
    if not np.all([np.all([np.all([np.all([wkeys[i][j][k][l] \
                                           == wkeys[0][j][0][0]\
                                           for l in range(nz)]) \
                                   for k in range(nph)]) \
                           for j in range(ncols)])\
                   for i in range(nrows)]):
        print('wkeys:\n', wkeys)
        raise RuntimeError('weights differ in same column')
    if not np.all([np.all([np.all([np.all([np.isclose(zvals[0][0][0][zi], 
                                                      zvals[i][j][pi][zi],
                                                      atol=3e-3, rtol=3e-3)\
                                           for zi in range(nz)])\
                                   for pi in range(nph)])\
                           for j in range(ncols)])\
                   for i in range(nrows)]):
        print('zvals:\n', zvals)
        raise RuntimeError('redshift values between panels or weights')
    if not np.all([np.all([np.all([np.all([pkeys[0][0][pi][0] \
                                           == pkeys[i][j][pi][zi] \
                                           for zi in range(nz)])\
                                   for pi in range(nph)])\
                           for j in range(ncols)])\
                   for i in range(nrows)]):
        print('pkeys:\n', pkeys)
        raise RuntimeError('physics models different between panels')
    
    styleargs_lines = {'med': {'linewidth': 2., 'linestyle': 'solid'},
                       'sct': {'linewidth': 1., 'linestyle': 'dashed'},
                       'sct_med': {'elinewidth': 1., 'linestyle': 'none'},
                       }
    styleargs_shade = {'med': {'linewidth': 2., 'linestyle': 'solid'},
                       'sct': {'alpha': 0.2},
                       'sct_med': {'elinewidth': 1., 'linestyle': 'none'},
                      }

    for ri in range(nrows):
        for ci in range(ncols):
            ax = axes[ri][ci]
            ykey = ykeys[ri][ci][0][0] # already checked equal in panel, row
            zval = zvals[ri][ci][0] # already checked equal in panels
            wkey = wkeys[ri][ci][0][0] # checked same in columns
            _pkeys = pkeys[ri][ci]
            _rvirs = rvirs[ri][ci]
            _rvals = rvals[ri][ci]
            _yvals = yvals[ri][ci]
            _hists_raw = hists_raw[ri][ci]
            
            _yvals = [yv - np.log10(ynormopts[ykey]) for yv in _yvals]
            if ri == nrows - 1:
                ax.set_xlabel(xlabel, fontsize=fontsize)
                labelx = True
            else:
                labelx = False
            if ri == 0:
                ax.set_title(wlabopts[wkey])
            if ci == 0:
                ax.set_ylabel(ylabopts[ykey], fontsize=fontsize)
                labely = True
            else:
                labely = False
            ax.tick_params(which='both', direction='in', top=True,
                           left=True, bottom=True, right=True,
                           labelbottom=labelx, labelleft=labely,
                           labelsize=fontsize - 1.)
            ax.set_xscale('log')
            ax.grid(visible=True)

            pv = np.array([0.1, 0.5, 0.9])
            zlab = f'$z={min(zval):.1f} \\endash {max(zval):.1f}$'
            rvir_min = np.inf
            rvir_max = -np.inf
            rc_all = None
            for pi, (__rvir, __rv, __yv, __hist_raw, __pkey) in\
                    enumerate(zip(_rvirs, _rvals, _yvals, _hists_raw, 
                                  _pkeys)):
                plos = []
                phis = []
                pmeds = []
                rc_all = None
                for rvir, rv, yv, hist_raw, pkey in \
                        zip(__rvir, __rv, __yv, __hist_raw, __pkey):
                    rvir_min = min(rvir_min, rvir)
                    rvir_max = max(rvir_max, rvir)
                    plo, pmed, phi = pu.percentiles_from_histogram(
                        10**hist_raw, yv, axis=1, percentiles=pv)
                    plos.append(plo)
                    phis.append(phi)
                    pmeds.append(pmed)
                    rc = 0.5 * (rv[1:] + rv[:-1])
                    if rc_all is None:
                        rc_all = rc
                    elif not np.allclose(rc_all, rc):
                        msg = 'r values mismatched in same panel, weight'
                        raise RuntimeError(msg)
                plos = np.array(plos)
                phis = np.array(phis)
                pmeds = np.array(pmeds)
                
                medmed = np.median(pmeds, axis=0)
                maxmed = np.max(pmeds, axis=0)
                minmed = np.min(pmeds, axis=0)
                mscthi = np.median(phis, axis=0)
                msctlo = np.median(plos, axis=0)
                deltahi = maxmed - medmed
                deltalo = medmed - minmed

                pstyle = 'lines'
                color = pcolors[pkey]
                medlab = 'med. of ' + zlab + ' med., ' + pkey
                sctlab_lines = 'med. of ' + zlab + ' 10, 90%, ' + pkey
                sctlab_shade = 'med. of ' + zlab + ' 10-90%, ' + pkey
                sct_medlab = zlab + ' med. range, ' + pkey
                if pstyle == 'lines':
                    ax.plot(rc_all, medmed, color=color, label=medlab,
                            **styleargs_lines['med'])
                    ax.plot(rc_all, mscthi, color=color, label=sctlab_lines,
                            **styleargs_lines['sct'])
                    ax.plot(rc_all, msctlo, color=color, 
                            **styleargs_lines['sct'])
                    ax.errorbar(rc_all, medmed, color=color, label=sct_medlab,
                                errorevery=(pi, nph), xerr=None, 
                                yerr=(deltalo, deltahi),
                                **styleargs_lines['sct_med'])
                elif pstyle == 'shade':
                    ax.plot(rc, medmed, color=color, label=medlab,
                            **styleargs_shade['med'])
                    ax.fill_between(rc, mscthi, msctlo, color=color, 
                                    **styleargs_shade['sct'],
                                    label=sctlab_shade)
                    ax.errorbar(rc_all, medmed, color=color, label=sct_medlab,
                                errorevery=(pi, nph), xerr=None, 
                                yerr=(deltalo, deltahi),
                                **styleargs_shade['sct_med'])
            if markrvir:
                ax.axvspan(rvir_min, rvir_max, alpha=0.2, color='gray')
            if ykey[0] == 'sim-direct' and \
                    ykey[1].startswith('ElementAbundance'):
                ymin, ymax = ax.get_ylim()
                ymin = max(ymin, -10.)
                ax.set_ylim(ymin, ymax)

        ylims = [axes[ri][i].get_ylim() for i in range(ncols)]
        ymin = min([ylim[0] for ylim in ylims])
        ymax = max([ylim[1] for ylim in ylims])
        [axes[ri][i].set_ylim((ymin, ymax)) for i in range(ncols)]
    handles, _ = axes[0][0].get_legend_handles_labels()
    lax.set_axis_off()
    lax.legend(handles=handles, fontsize=fontsize - 1., 
               ncol=int(np.floor(0.5 * ncols)),
               loc='upper center', bbox_to_anchor=(0.5, 0.65))

    if figtitle is not None:
        fig.suptitle(figtitle, fontsize=fontsize)

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotset3dprof_physcomp_weightvalgrid(fileset):
    if fileset == 'clean_set1-2':
        filedir = '/Users/nastasha/ciera/profiles/fire/clean_set1_set2/'
        filen_template = ('hist_{yq}_r3D_by_{wq}_{simname}'
                          '_snap{snap}_bins1_v1.hdf5')
        simnames = [[('m13h113_m3e5_MHDCRspec1_fire3_fireBH'
                      '_fireCR1_Oct252021_crdiffc1_sdp1e-4'
                      '_gacc31_fa0.5_fcr1e-3_vw3000'),
                      ('m13h113_m3e4_MHD_fire3_fireBH'
                       '_Sep182021_hr_crdiffc690_sdp1e-4'
                       '_gacc31_fa0.5'),
                      ('m13h113_m3e5_MHD_fire3_fireBH'
                       '_Sep182021_crdiffc690_sdp1e10'
                       '_gacc31_fa0.5'),
                     ],
                     [('m13h206_m3e5_MHDCRspec1_fire3_fireBH'
                       '_fireCR1_Oct252021_crdiffc1_sdp1e-4'
                       '_gacc31_fa0.5_fcr1e-3_vw3000'),
                      ('m13h206_m3e4_MHD_fire3_fireBH'
                       '_Sep182021_hr_crdiffc690_sdp3e-4'
                       '_gacc31_fa0.5'),
                      ('m13h206_m3e5_MHD_fire3_fireBH'
                       '_Sep182021_crdiffc690_sdp1e10'
                       '_gacc31_fa0.5'),
                     ],
                     [('m12f_m6e4_MHDCRspec1_fire3_fireBH'
                       '_fireCR1_Oct252021_crdiffc1_sdp1e-4'
                       '_gacc31_fa0.5_fcr1e-3_vw3000'),
                      ('m12f_m7e3_MHD_fire3_fireBH'
                       '_Sep182021_hr_crdiffc690_sdp2e-4'
                       '_gacc31_fa0.5'),
                      ('m12f_m7e3_MHD_fire3_fireBH'
                       '_Sep182021_hr_crdiffc690_sdp1e10'
                       '_gacc31_fa0.5'),
                     ],
                    ] 
        simlabs = ['m13h113', 'm13h206', 'm12f']

        srst = ['m13h113_m3e5', 'm13h206_m3e5', 'm12f_m6e4']
        hrst = ['m13h113_m3e4', 'm13h206_m3e4', 'm12f_m7e3']
        snaps_sr = [50, 49, 48, 47, 46, 45]
        snaps_hr = [258, 240, 224, 210, 197, 186]
        #redshifts = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        simlen = len(simnames)
        fillkw_valss = [[{'yq': 'Density'},
                         {'yq': 'Temperature'},
                         {'yq': 'Neon'},
                         {'yq': 'Oxygen'},
                         {'yq': 'Magnesium'},
                        ],
                        [{'yq': 'Density'},
                         {'yq': 'Temperature'},
                         {'yq': 'Neon'},
                        ],
                        [{'yq': 'Density'},
                         {'yq': 'Temperature'},
                         {'yq': 'Oxygen'},
                        ],
                        [{'yq': 'Density'},
                         {'yq': 'Temperature'},
                         {'yq': 'Magnesium'},
                        ]
                       ] * simlen
        fillkw_weightss = [[{'wq': 'Mass'}, 
                            {'wq': 'Volume'}, 
                            {'wq': 'H1'},
                           ],
                           [{'wq': 'Mass'}, 
                            {'wq': 'Volume'}, 
                            {'wq': 'Ne8'},
                           ],
                           [{'wq': 'Mass'}, 
                            {'wq': 'Volume'}, 
                            {'wq': 'O6'},
                           ],
                           [{'wq': 'Mass'}, 
                            {'wq': 'Volume'}, 
                            {'wq': 'Mg10'},
                           ],
                          ] * simlen
        weightlabs = ['M_V_HI', 'M_V_Ne8', 'M_V_O6', 'M_V_Mg10']
        fillkw_zs = [[[{'snap': snhr} \
                            if np.any([sim.startswith(st) for st in hrst])\
                        else {'snap': snsr} \
                            if np.any([sim.startswith(st) for st in srst])\
                        else None \
                       for snhr, snsr in zip(snaps_hr, snaps_sr)] \
                      for sim in simlist] \
                     for simlist in simnames for i in weightlabs]
        filen_templates = \
            [filedir + filen_template \
             for i in range(simlen) for weight in weightlabs]
        outdir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/3dprof/'
        outtemplate = 'prof3D_physcomp_{wv}_{simlab}.pdf'
        outnames = [outdir + outtemplate.format(
                        simlab=simlab.replace(' ', '_'),
                        wv=weight)\
                    for simlab in simlabs for weight in weightlabs]
        figtitles = [simlab for simlab in simlabs for weight in weightlabs]
        fillkw_physs = [[{'simname': sim} for sim in simlist]\
                         for simlist in simnames for i in weightlabs]

    for filen_template, fillkw_weights, fillkw_phys, fillkw_vals, fillkw_z, \
            figtitle, outname in zip(filen_templates, fillkw_weightss, 
                                     fillkw_physs, fillkw_valss, fillkw_zs, 
                                     figtitles, outnames):
        plot3dprof_physcomp_weightvalgrid(filen_template, fillkw_weights, 
                                          fillkw_vals, fillkw_z, fillkw_phys,
                                          figtitle=figtitle, outname=outname)
        
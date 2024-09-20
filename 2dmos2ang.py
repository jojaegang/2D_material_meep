# %%
import meep as mp
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy
# from meep.materials import Ag
def def_MoS2(MoS2_thick, eps_background):

    # default unit length is 1 μm
    um_scale = 1.0
    pixel_thickness = MoS2_thick
    # pixel_thickness = 0.61/1000

    monolayer_thickness = 0.61/1000
    # conversion factor for eV to 1/μm [=1/hc]
    eV_um_scale = um_scale / 1.23984193

    # ------------------------------------------------------------------
    # MoS2 monolayer from https://doi.org/10.1002/adpr.202000180
    # wavelength range: 0.38 - 0.66 μm

    # MoS2_freq_range = mp.FreqRange(min=um_scale/0.7, max=um_scale / 0.38)

    MoS2_eps_infty = 4.83291507293635
    MoS2_frq1 = 1.64797111642578
    MoS2_gamma1 = 0.108826471219756
    MoS2_sig1 = 0.417364996349752
    MoS2_frq2 = 1.53445492221022
    MoS2_gamma2 = 0.0380034929217783
    MoS2_sig2 = 0.126706036520689
    MoS2_frq3 = 2.32438579875489
    MoS2_gamma3 = 0.264618881073392
    MoS2_sig3 = 1.38265645177994
    MoS2_frq4 = 2.54804525191946
    MoS2_gamma4 = 0.822171117909694
    MoS2_sig4 = 3.10838595757667

    MoS2_eps_infty = 1 + (monolayer_thickness/pixel_thickness)*(MoS2_eps_infty-1)
    MoS2_sig1 = MoS2_sig1*(monolayer_thickness/pixel_thickness)
    MoS2_sig2 = MoS2_sig2*(monolayer_thickness/pixel_thickness)
    MoS2_sig3 = MoS2_sig3*(monolayer_thickness/pixel_thickness)
    MoS2_sig4 = MoS2_sig4*(monolayer_thickness/pixel_thickness)

    MoS2_susc = [
        mp.LorentzianSusceptibility(frequency=MoS2_frq1, gamma=MoS2_gamma1, sigma_diag=MoS2_sig1*mp.Vector3(1,0,1)),
        mp.LorentzianSusceptibility(frequency=MoS2_frq2, gamma=MoS2_gamma2, sigma_diag=MoS2_sig2*mp.Vector3(1,0,1)),
        mp.LorentzianSusceptibility(frequency=MoS2_frq3, gamma=MoS2_gamma3, sigma_diag=MoS2_sig3*mp.Vector3(1,0,1)),
        mp.LorentzianSusceptibility(frequency=MoS2_frq4, gamma=MoS2_gamma4, sigma_diag=MoS2_sig4*mp.Vector3(1,0,1)),
    ]

    MoS2 = mp.Medium(epsilon_diag=mp.Vector3(MoS2_eps_infty+(eps_background-1),eps_background,MoS2_eps_infty+(eps_background-1)), 
                     E_susceptibilities=MoS2_susc, 
                    #  valid_freq_range=MoS2_freq_range,
                     )
    return MoS2


# %%
def flux_sim(freq, resolution, MoS2, MoS2_thick, source_component, ang, eps_background):
#center position half pixel, or pixel???, thickness one pixel or two??? comparison with TMM
    # when the # of pixel is odd, integer center. when the # of pixel is even, half integer.

    Lx = 0.1 # fix
    Ly = 3
    Lpml = 0.5
    cell_size = mp.Vector3(Lx, Ly)
    dimensions = 2
    pml_layers = [mp.PML(Lpml, direction=mp.Y)]
    fwidth = 0.04

    kx_lamda = 2*np.pi*np.sin(ang*np.pi/180)*freq*abs(np.sqrt(eps_background))

    src = mp.GaussianSource(frequency=freq, fwidth=fwidth, is_integrated=True)
    sources = []

    for i in range(round(resolution*Lx)):
        # if i%3==1 or i%3==2:
        #     continue
        sources.append(mp.Source(src,
                        component=source_component,
                        center=mp.Vector3(x = Lx*(i/(Lx*resolution)-0.5)+0.5/resolution, y = -0.5*Ly+Lpml+0.1),
                        # size = mp.Vector3(),
                        amplitude = (np.exp(1j*kx_lamda*Lx*(i/(Lx*resolution)-0.5))) *1/(resolution*Lx),
                        ))

    sim = mp.Simulation(
        cell_size=cell_size,
        sources=sources,
        resolution=resolution,
        boundary_layers=pml_layers,
        default_material = mp.Medium(epsilon = eps_background),
        # geometry=geometry,
        dimensions=dimensions,
        force_complex_fields=False,
        k_point=mp.Vector3(kx_lamda/2/np.pi,0,0),
        eps_averaging=False,
    )
    sim_time = 10 # fix

    l_monitor_pos = -0.5*Ly+Lpml+0.5
    r_monitor_pos = +0.5*Ly-Lpml-0.5
    print(l_monitor_pos)
    print(r_monitor_pos)


    flux_l = mp.FluxRegion(
        center=mp.Vector3(y=l_monitor_pos), 

        size=mp.Vector3(x = Lx, y = 0)
        )
    n_freq  = 1
    flux_f_width = 0.00
    flux_l = sim.add_flux(freq,flux_f_width,n_freq,flux_l,)


    # dft_fields = sim.add_dft_fields([source_component],
    #                                 freq,0,1,
    #                                 center=mp.Vector3(),
    #                                 size=mp.Vector3(Lx, Ly),
    #                                 yee_grid=True)

    sim.run(until_after_sources= sim_time)

    freq_array = mp.get_flux_freqs(flux_l)
    empty_flux_l_array = mp.get_fluxes(flux_l)
    empty_flux_l_data = sim.get_flux_data(flux_l)

    geometry = [
        mp.Block(mp.Vector3(x = Lx, y=MoS2_thick), center = mp.Vector3(x=0, y=0./resolution), material = MoS2),
                ]

    sim.reset_meep()

    sim = mp.Simulation(
        cell_size=cell_size,
        sources=sources,
        resolution=resolution,
        boundary_layers=pml_layers,
        default_material = mp.Medium(epsilon = eps_background),
        geometry=geometry,
        dimensions=dimensions,
        force_complex_fields=False,
        k_point=mp.Vector3(kx_lamda/2/np.pi,0,0),
        eps_averaging=False,
    )
    flux_l = mp.FluxRegion(
        center=mp.Vector3(y=l_monitor_pos), 
        size=mp.Vector3(x = Lx, y = 0)
        )
    flux_r = mp.FluxRegion(
        center=mp.Vector3(y=r_monitor_pos), 
        size=mp.Vector3(x = Lx, y = 0)

        # size=mp.Vector3()
        )

    print(l_monitor_pos)
    print( r_monitor_pos)
    flux_r = sim.add_flux(freq,flux_f_width,n_freq,flux_r,)
    flux_l = sim.add_flux(freq,flux_f_width,n_freq,flux_l,)
    # flux_l_data = sim.get_flux_data(flux_l)
    sim.load_minus_flux_data(flux_l, empty_flux_l_data)

    sim.run(until_after_sources= sim_time)


    freq_array = mp.get_flux_freqs(flux_r)
    flux_r_array = mp.get_fluxes(flux_r)
    flux_l_array = mp.get_fluxes(flux_l)
    return [freq_array, flux_l_array, flux_r_array, empty_flux_l_array]

# %%
resolution = 50 #fix in all the angles
MoS2_thick = 1.0/resolution
# eps_background = 1
eps_background = 1.4**2

MoS2 = def_MoS2(MoS2_thick,eps_background)  #---------------------input



component = 'Ez'
source_component = mp.Ez # Ez, Hz ---------------------input

init_wavelength = 0.38
final_wavelength = 0.68

for ang in np.array([20, 40, 60, 80]):
    for freq in 1./np.linspace(init_wavelength, final_wavelength, 30):
        [freq_array, flux_l_array, flux_r_array, empty_flux_l_array] = flux_sim(freq, resolution, MoS2, MoS2_thick, source_component, ang, eps_background)
        trans_data = np.vstack((np.ones(np.size(freq_array))/freq_array, np.array(flux_r_array)/np.array(empty_flux_l_array)))
        trans_data = np.transpose(trans_data)
        refl_data = np.vstack((np.ones(np.size(freq_array))/freq_array, -np.array(flux_l_array)/np.array(empty_flux_l_array)))
        refl_data = np.transpose(refl_data)
        abs_data = np.vstack((np.ones(np.size(freq_array))/freq_array, 1-np.array(flux_r_array)/np.array(empty_flux_l_array)+np.array(flux_l_array)/np.array(empty_flux_l_array)))
        abs_data = np.transpose(abs_data)

        if freq == 1/init_wavelength:
            ang_trans_data  = trans_data
            ang_abs_data  = abs_data
            ang_refl_data  = refl_data
        ang_trans_data = np.concatenate((ang_trans_data,trans_data),axis = 0)
        ang_abs_data = np.concatenate((ang_abs_data,abs_data),axis = 0)
        ang_refl_data = np.concatenate((ang_refl_data,refl_data),axis = 0)
        print("ang: "+str(ang)+", wavelength: " + str(round(1000/freq)))

    np.savetxt('MoS2_trans_'+str(component)+'_'+str(ang)+'.txt', ang_trans_data, encoding='utf8', delimiter=' ' )
    np.savetxt('MoS2_abs_'+str(component)+'_'+str(ang)+'.txt', ang_abs_data, encoding='utf8', delimiter=' ' )
    np.savetxt('MoS2_refl_'+str(component)+'_'+str(ang)+'.txt', ang_refl_data, encoding='utf8', delimiter=' ' )

    

# %%
# plt.plot(1000*np.ones(n_freq)/freq_array, np.array(flux_r_array)/np.array(empty_flux_l_array))
# plt.xlabel("Wavelength (nm)")
# plt.ylabel("Transmittance")
# plt.show()
# data = np.vstack((np.ones(n_freq)/freq_array, np.array(flux_r_array)/np.array(empty_flux_l_array)))
# data = np.transpose(data)
# np.savetxt('1D_MoS2_trans_arr.txt', data, encoding='utf8', delimiter=', ' )

# plt.plot(1000*np.ones(n_freq)/freq_array,  1-np.array(flux_r_array)/np.array(empty_flux_l_array)+np.array(flux_l_array)/np.array(empty_flux_l_array))
# # plt.plot(np.ones(n_freq)/freq_array, -np.array(flux_l_array))
# plt.ylim([0,0.2])
# plt.xlabel("Wavelength (nm)")
# plt.ylabel("Absorption")
# plt.show()
# # plt.plot(np.ones(n_freq)/freq_array, (np.array(flux_l_array)+np.array(free_flux_r_array))/np.array(free_flux_r_array))
# sim.get_flux_data(flux_l)
# plt.plot(1000*np.ones(n_freq)/freq_array, -np.array(flux_l_array)/np.array(empty_flux_l_array))
# # plt.plot(np.ones(n_freq)/freq_array, -np.array(flux_l_array))


# plt.xlabel("Wavelength (nm)")
# plt.ylabel("Refelectance")
# plt.show()

# data = np.vstack((np.ones(n_freq)/freq_array, 1-np.array(flux_r_array)/np.array(empty_flux_l_array)+np.array(flux_l_array)/np.array(empty_flux_l_array)))
# data = np.transpose(data)
# np.savetxt('1D_MoS2_absorb_arr.txt', data, encoding='utf8', delimiter=', ' )

# data = np.vstack((np.ones(n_freq)/freq_array, -np.array(flux_l_array)/np.array(empty_flux_l_array)))
# data = np.transpose(data)
# np.savetxt('1D_MoS2_reflect_arr.txt', data, encoding='utf8', delimiter=', ' )



def main():

    import pygmo as pg
    import pykep as pk
    
    # Use the JPL low precision ephemerides
    from pykep.planet import jpl_lp 
    
    from matplotlib import pyplot as plt
    import numpy as np
    from scipy import array

    #Pykep makes use of lambert transfers to study interplanetary trajectories by producing Pork-Chop plots

    #Start by sampling departure epochs and transfer times every 15 days and solve the lambert problem in a large defined grid
    
    start_epochs = np.arange(4000.0,10000.0,15.0)
    duration = np.arange(180.0,470.0,15.0)
    
    # Define the start and end planets of the inteplanetary trajectory
    start_planet = jpl_lp('earth')
    end_planet = jpl_lp('mars')

    data = list()
    for start in start_epochs:
        row = list()
        for T in duration:
            r1,v1 = start_planet.eph(pk.epoch(start))
            r2,v2 = end_planet.eph(pk.epoch(start+T))
            l = pk.lambert_problem(r1,r2,T*60*60*24, start_planet.mu_central_body)
            DV1 = np.linalg.norm(array(v1)-array(l.get_v1()[0]))
            DV2 = np.linalg.norm(array(v2)-array(l.get_v2()[0]))
            DV1 = max([0,DV1-4000])
            DV = DV1+DV2
            row.append(DV)
        data.append(row)

    #Extract best sol and relative epochs
    minrows = [min(l) for l in data]
    i_idx = np.argmin(minrows)
    j_idx = np.argmin(data[i_idx])
    best = data[i_idx][j_idx]
    print('Best DV: ',best)
    print('Launch epoch (MJD2000): ',start_epochs[i_idx])
    print('Duration (days): ',duration[j_idx])
    duration_pl, start_epochs_pl = np.meshgrid(duration, start_epochs)
    CP1 = plt.contourf(start_epochs_pl,duration_pl,array(data),levels=list(np.linspace(best,5000,10)))
    plt.colorbar(CP1).set_label('△V km/s')   
    
    plt.title(f'{start_planet.name} - {end_planet.name} Total △V Requirements'.title())
    plt.xlabel('Launch Date (MJD2000)')
    plt.ylabel('Mission Duration (days)')
    plt.show()


    #ZOOM IN TO BEST EPOCH - sample with finer resolution
    start_epochs = np.arange(8950.0,9100.0,1.0)
    duration = np.arange(200.0,450.0,1.0)
    data=list()
    for start in start_epochs:
        row = list()
        for T in duration:
            r1,v1 = start_planet.eph(pk.epoch(start))
            r2,v2 = end_planet.eph(pk.epoch(start+T))
            l = pk.lambert_problem(r1,r2,T*60*60*24, start_planet.mu_central_body)
            DV1 = np.linalg.norm(array(v1)-array(l.get_v1()[0]))
            DV2 = np.linalg.norm(array(v2)-array(l.get_v2()[0]))
            DV1 = max([0,DV1-4000])
            DV = DV1+DV2
            row.append(DV)
        data.append(row)

    minrows = [min(l) for l in data]
    i_idx = np.argmin(minrows)
    j_idx = np.argmin(data[i_idx])
    best = data[i_idx][j_idx]
    print('Best DV: ',best)
    print('Launch epoch (MJD2000): ',start_epochs[i_idx])
    print('Duration (days): ',duration[j_idx])
    duration_pl2, start_epochs_pl2 = np.meshgrid(duration, start_epochs)

    CP2 = plt.contourf(start_epochs_pl2,duration_pl2,array(data),levels=list(np.linspace(best,5000,10)))
    plt.colorbar(CP2).set_label('△V km/s')  
    plt.title(f'{start_planet.name} - {end_planet.name} Total △V Requirements'.title())
    plt.xlabel('Launch Date (MJD2000)')
    plt.ylabel('Mission Duration (days)')
    plt.show()


if __name__ == "__main__":
	main()

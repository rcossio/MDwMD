# This script is to calculate the diffusion coefficient from the covariance base estimator (CVE)
# as indicated by Bullerjahn et al. 2020 (doi: 10.1063/5.0008312)

import numpy as np
import sys
import matplotlib.pyplot as plt

def load_data(filepath):
    """Load the COM data from an xvg file."""
    try:
        data = np.loadtxt(filepath, comments=["@", "#"])
        time = data[:, 0]
        com = data[:, 1:4]
        return time, com
    except Exception as e:
        print(f"Failed to load data: {e}")
        sys.exit(1)


def main():
    filepath = 'prod_com.xvg'
    time, com = load_data(filepath)
    total_frames = com.shape[0]
    delta_t = time[1] - time[0]

    print(f"Reading {total_frames} frames with a time step of {delta_t} ps")

    analysed_delta_t, estimated_D, estimated_var_D = [], [], []


    for i in range(1,1000+1):
        sub_delta_t = i*delta_t
        x = com[::i]  #every frame
        dx = x[1:] - x[:-1]  #displacement data frame

        # Calculate a2
        a2 = dx[1:]*dx[:-1]  #covariance data frame
        a2 = -2* np.mean(a2, axis=0)     #calc a2 per dimension x,y,z
        a2_3D = np.sum(a2)  #sum of all dimensions

        # Calculate s2
        s2=dx**2    #variance data frame
        s2=np.mean(s2, axis=0)-a2  #calc s2 per dimension x,y,z
        s2_3D=np.sum(s2)  #sum of all dimensions

        # Calculate D
        D=1e6*s2_3D/6/sub_delta_t  #diffusion coefficient in um^2/s

        # Calculate variances
        N = len(dx)
        var_a2  = (7*a2**2+8*a2*s2+4*s2**2)/(N-1)
        var_a2 -=  (2*a2**2)/(N-1)**2

        var_s2  = 4*(a2*s2+s2**2)/(N-1)
        var_s2 += 2*(a2**2+s2**2)/N
        var_s2 += (5*a2**2+4*a2*s2)/N/(N-1)
        var_s2 -= a2**2/(N-1)**2
        var_s2 -= a2**2/N**2/(N-1)**2

        var_D=np.sum(var_s2)/(6*sub_delta_t)**2
        var_D=1e6*np.sqrt(var_D)

        #if (a2_3D>0) or ((s2_3D/a2_3D)>1.0):
        analysed_delta_t.append(sub_delta_t)
        estimated_D.append(D)
        estimated_var_D.append(var_D)
        print(f"Timestep: {sub_delta_t:.0f} ps, N={len(x)}, D={D:.1f}({var_D:.1f}), a2={a2_3D:.5f}, s2={s2_3D:.5f}, s2/a2={s2_3D/a2_3D:.5f}")
        
    #plotting
    plt.figure(figsize=(10, 6))
    plt.plot(analysed_delta_t, estimated_D, label='Diffusion Coefficient')
    plt.fill_between(analysed_delta_t, np.array(estimated_D)-np.array(estimated_var_D), np.array(estimated_D)+np.array(estimated_var_D), alpha=0.3)
    plt.title('D(t)')
    plt.grid(True)
    plt.xlabel('Time (ps)')
    plt.ylabel('DiffusionMSD (um$^{2}$/s)')
    plt.legend()
    plt.savefig('diffusion_CVE.png')



if __name__ == "__main__":
    main()

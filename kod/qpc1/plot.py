# pyright: reportUndefinedVariable=false
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({"text.usetex":True, "font.size": 16})

def wczytaj(*nazwy): # funkcja do wczytywania plikow od clauda
    for nazwa in nazwy:
        globals()[nazwa] = np.loadtxt(nazwa + ".csv", delimiter=',')

# wczytanie plików
wczytaj("dyspfile","kxfile","lfile","ufile","poprzecznefile")

Ny = len(lfile)

uRe = ufile[:,0].reshape(Ny,Ny)
uIm = ufile[:,1].reshape(Ny,Ny)
poprzecznefile = poprzecznefile.astype(bool)

size = 8

print("Minima pasm:")
for j in range(len(dyspfile[0])):
    if j<4: print(dyspfile[:,j].min())

# relacja dyspersji
plt.figure(figsize=(1.25*size,size))
for j in range(len(dyspfile[0])):
    plt.plot(kxfile,dyspfile[:,j],c="tab:blue")
plt.grid(ls=":")
plt.xlabel(rf"$k_x\ \left[\frac{{1}}{{nm}}\right]$")
plt.ylabel(rf"$E\ \left[eV\right]$")
plt.tight_layout()
plt.savefig("dyspfile.png",dpi=150)
plt.close()

# wszystkie wartosci wlasne upw
plt.figure(figsize=(size,size))
plt.scatter(range(len(lfile)),lfile)
plt.axhline(1.0,ls='--',c='k')
plt.savefig("lfile.png")
plt.close()

# wszystkie wektory wlasne upw
# Re
plt.figure(figsize=(1.25*size,size))
plt.imshow(uRe,origin='lower',cmap='bwr')
plt.colorbar()
plt.tight_layout()
plt.savefig("uRe.png")
plt.close()
# Im
plt.figure(figsize=(1.25*size,size))
plt.imshow(uIm,origin='lower',cmap='bwr')
plt.colorbar()
plt.tight_layout()
plt.savefig("uIm.png")
plt.close()

# mody poprzeczne
# Re
licznik=0
plt.figure(figsize=(size,size))
for j in range(Ny):
    if (poprzecznefile[j]): 
        plt.plot(uRe[j],label=rf"$\Re(u_{{-,{licznik}}})$")
        plt.plot(uIm[j],label=rf"$\Im(u_{{-,{licznik}}})$")
        licznik+=1
plt.legend()
plt.savefig("popU.png")
plt.close()

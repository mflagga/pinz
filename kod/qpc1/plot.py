# pyright: reportUndefinedVariable=false
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 20,
    "text.latex.preamble": r"""
        \usepackage{siunitx}
        \DeclareSIUnit{\arbitrary}{a.u.}
    """,
})

def wczytaj(*nazwy): # funkcja do wczytywania plikow od clauda
    for nazwa in nazwy:
        globals()[nazwa] = np.loadtxt(nazwa + ".csv", delimiter=',')

# wczytanie plików
wczytaj("dyspfile","kxfile","lfile","ufile","poprzecznefile","kx2vfile","misc","psifile")

E=misc[0]
Ny = int(misc[1])
Nx = int(misc[2])
xmin=misc[3]
ymin=misc[4]
liczba=int(misc[5])
kx2vfile = np.atleast_2d(kx2vfile)
psifile = np.atleast_2d(psifile)
uRe = ufile[:,0].reshape(Ny,Ny)
uIm = ufile[:,1].reshape(Ny,Ny)
poprzecznefile = poprzecznefile.astype(bool)
Evec=np.ones(len(kx2vfile[:,0]))*E
dx=-2*xmin/(Nx+1)
dy=-2*ymin/(Ny+1)
# for l in range(liczba):
#     psifile[l].reshape(Nx,Ny)

x=np.linspace(xmin+dx,-xmin-dx,Nx)
y=np.linspace(ymin+dy,-ymin-dy,Ny)

size = 8

# print("Minima pasm:")
# for j in range(len(dyspfile[0])):
#     if j<4: print(dyspfile[:,j].min())

# relacja dyspersji
plt.figure(figsize=(1.25*size,size))
for j in range(len(dyspfile[0])):
    plt.plot(kxfile,dyspfile[:,j],c="tab:blue")
plt.scatter(kx2vfile[:,0],Evec,c='r')
plt.scatter(-kx2vfile[:,0],Evec,c='r')
plt.axhline(E,ls='--',c='k',label=rf"$E=\qty{{{E}}}{{\electronvolt}}$")
plt.grid(ls=":")
plt.legend()
plt.xlabel(r"$k_x\ [\unit{\per\nano\meter}]$")
plt.ylabel(r"$E\ [\unit{\electronvolt}]$")
plt.tight_layout()
plt.savefig("dyspfile.png",dpi=150)
plt.close()

# wszystkie wartosci wlasne upw
plt.figure(figsize=(size,size))
plt.scatter(range(len(lfile)),lfile)
plt.axhline(1.0,ls='--',c='k')
plt.xlabel(r"$n$")
plt.ylabel(r"$|\lambda_n^-|$")
plt.tight_layout()
plt.savefig("lfile.png",dpi=150)
plt.close()

# wszystkie wektory wlasne upw
# Re
Remax = np.abs(uRe).max()
plt.figure(figsize=(1.25*size,size))
plt.imshow(uRe,origin='lower',cmap='bwr',vmin=-Remax,vmax=Remax)
plt.colorbar()
plt.title(r"$\Re(u^-)\ [\unit{\nano\meter\tothe{-\frac{1}{2}}}]$")
plt.tight_layout()
plt.savefig("uRe.png",dpi=150)
plt.close()
# Im
Immax = np.abs(uIm).max()
plt.figure(figsize=(1.25*size,size))
plt.imshow(uIm,origin='lower',cmap='bwr',vmin=-Immax,vmax=Immax)
plt.colorbar()
plt.title(r"$\Im(u^-)\ [\unit{\nano\meter\tothe{-\frac{1}{2}}}]$")
plt.tight_layout()
plt.savefig("uIm.png",dpi=150)
plt.close()

# mody poprzeczne
# Re
licznik=0
plt.figure(figsize=(size,size))
for j in range(Ny):
    if (poprzecznefile[j]): 
        plt.plot(y,uRe[j],label=rf"$\Re(u_{{-,{licznik}}})$")
        plt.plot(y,uIm[j],label=rf"$\Im(u_{{-,{licznik}}})$")
        licznik+=1
plt.xlabel(r"$y\ [\unit{\nano\meter}]$")
plt.ylabel(r"$u_{-,n}\ [\unit{\nano\meter\tothe{-\frac{1}{2}}}]$")
plt.legend()
plt.tight_layout()
plt.savefig("popU.png",dpi=150)
plt.close()

plt.figure(figsize=(size,size))
plt.imshow(psifile[1,:].reshape(Nx,Ny).T,origin='lower',extent=[xmin, -xmin, ymin, -ymin])
plt.colorbar()
plt.savefig("psifile.png",dpi=150)

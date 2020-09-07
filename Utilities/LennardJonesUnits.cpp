#include <math.h>
#include "LennardJonesUnits.h"

void CLJUnits::defineUnits()
{
    // Divide real units with these to obtain value in reduced units (i.e. Lennard-Jones like units).
    // Units into and out of a program using LJ like units are in number of units of for example
    // m_dTimeUnit or m_dVelocity units.
    massUnit_ = 1.0;                      // [g/mol]
    distanceUnit_ = 1.0;                  // [AA]
    energyUnit_ = 1.0;                    // [kJ/mol]
    timeUnit_ = 100.0;                    // => [fs]
    velocityUnit_ = 0.01;                 // => [AA/fs]
    forceUnit_ = 1.0E10;                  // => [kN/mol]
    torqueUnit_ = energyUnit_;            // => [kJ/mol]
    tempUnit_ = 120.272;                  // => [K]
    pressUnit_ = 16388.2461;              // => [atm]
    chargeUnit_ = 0.0268283;              // => [units of |e|]
    volumeUnit_ = 1.0;                    // => [AA^3]
    buckinghamC_ = 1.0;                   // => [AA^6kJ/mol]
    harmonicBondK_ = 1.0;                 // => [(kJ/mol)/AA^2]
    
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Theory / Derivation (in LaTeX code)
    // -------------------------------------------
    // Let $\mathcal{M}$, $\mathcal{D}$ and $\mathcal{E}$ be the units of mass, distance and energy,
    // respectively. For a quantity $A$ in SI units, we denote the corresponding quantity in units of
    // $M$, $D$ and $E$ by $\tilde{A}$. Unit conversion factors are caligraphic or greek. Then, since
    // we require that
    // $$\tilde{E}=\tilde{m}\frac{d^2\tilde{x}}{d\tilde{t}^2}\Delta{\tilde{r}},$$
    // for energy $E$, mass $m$, distances $x$ and $\Delta r$, we have that
    // $$\frac{E}{\mathcal{E}}=\frac{m}{\mathcal{M}}\frac{d^2x/\mathcal{D}}{dt^2/\tau^2}\Delta{r}/\mathcal{D}.$$
    // Furthermore, since we also have that
    // $$E=m\frac{d^2x}{dt^2}\Delta{r},$$
    // we must have a derived unit of time
    // $$\tau = \sqrt{\frac{\mathcal{M}\mathcal{D}^2}{\mathcal{E}}}.$$
    // We can find the relationship for other derived units in a similar manner. For velocity
    // $$\frac{v}{\mathcal{V}}=\tilde{v}=\frac{dx/\mathcal{D}}{dt/\tau} \Rightarrow \mathcal{V}=\frac{\mathcal{D}}{\tau}.$$
    // For force
    // $$\tilde{E}=\frac{E}{\mathcal{D}}=\frac{F}{\mathcal{F}}\frac{\Delta r}{\mathcal{D}} \Rightarrow \mathcal{F}=\frac{\mathcal{E}}{\mathcal{D}}.$$
    // Torque has units of energy. The temperature unit is in units of energy, $\mathcal{E}$, s.t. internally
    // in the program using $\mathcal{M}$, $\mathcal{D}$, $\mathcal{E}$ units we have that the equipartition
    // principle reads $\tilde{E}_k=(1/2)\tilde{T}$. Thus,
    // $$\tilde{E}_k=\frac{1}{2}k_B T,$$
    // where $k_B$ is the Boltzmann constant. Therefore, if we want [$T$]=K, then [$k_B$]=(kJ/mol)/K
    // to make sure that [$\tilde{E}$]=kJ/mol. To achive this,
    // $$k_B=1.38064852\cdot 10^{-23}\cdot 10^{-3}\cdot N_A\text{kJ/K} = 0.00831446 \text{(kJ/mol)/K},$$
    // where $N_A$ is Avogadro's constant. We therefore have that $\tilde{T}=T/\mathcal{T}$, where
    // $$\mathcal{T}=(1/k_B)=120.272\text{K/(kJ/mol)}.$$
    // For pressure we have
    // $$\tilde{p}=\frac{p}{\mathcal{P}}=\frac{\tilde{F}}{\tilde{A}}=\frac{\tilde{F}}{\tilde{x}^2}=\frac{F/\mathcal{F}}{x^2/\mathcal{D}^2} \Rightarrow \mathcal{P}=\frac{\mathcal{F}}{\mathcal{D}^2}=\frac{\mathcal{E}}{\mathcal{D}^3}.$$
    // Charge unit is defined as
    // $$\mathcal{Q}=\sqrt{4\pi\epsilon_0\mathcal{D}\mathcal{E}},$$
    // such that $\tilde{Q}=Q/\mathcal{Q}$. For volume we have
    // $$\tilde{V}=\frac{V}{\nu}=\tilde{x}^3=\frac{x^3}{\mathcal{D}^3} \Rightarrow \nu=\mathcal{D}^3.$$
    // For Buckingham $C$ (often seen in units of \AA$^6$kJ/mol) we have
    // $$\tilde{C}=\frac{C}{\mathcal{C}}=\tilde{E}\tilde{x}^6=\frac{E}{\mathcal{E}}\frac{x^6}{\mathcal{D}^6} \Rightarrow \mathcal{C}=\mathcal{E}\mathcal{D}^6.$$
    //
    // We now want to calculate the number of fs in 1$\tau$ if $\mathcal{M}=1$g/mol, $\mathcal{D}=1$\AA~and $\mathcal{E}=1$kJ/mol:
    // $$1\tau=10^{-13}\sqrt{\frac{\text{kgm}^2}{\text{J}}}=10^{-13}\sqrt{\frac{\text{kgm}^2}{\text{kg}(\text{m}/\text{s}^2)\text{m}}}=10^{-13}\text{s}=10^2\text{fs},$$
    // thus yielding 1$\tau$=$10^2$fs, giving 1fs=$10^{-2}\tau$. Therefore, to give a value of 3fs
    // to a program in our LJ like units we must deliver a value of $3\cdot 10^{-2}$. Similarly,
    // $$1\mathcal{V}=\frac{\mathcal{D}}{\tau}=\frac{1\text{\AA}}{10^2\text{fs}}=10^{-2}\text{\AA/fs}\Rightarrow 1\text{\AA/fs}=10^2\mathcal{V},$$
    // $$1\mathcal{F}=\frac{\mathcal{E}}{\mathcal{D}}=\frac{1\text{kJ/mol}}{1\text{\AA}}=10^{10}\text{kN/mol} \Rightarrow 1\text{kN/mol}=10^{-10}\mathcal{F},$$
    // $$1\mathcal{P}=\frac{\mathcal{E}}{\mathcal{D}^3}=\frac{\text{kJ/mol}}{\text{\AA}^3}=10^{33}\frac{\text{N/mol}}{\text{m}^2}=\frac{10^{33}}{6.022\dots 10^{23}}\text{Pa}=16388.2\text{atm} \Rightarrow 1\text{atm}=\frac{1}{16388.2}\mathcal{P},$$
    // \begin{align*} 1\mathcal{Q}&=\sqrt{4\pi\epsilon_0\mathcal{D}\mathcal{E}}=\sqrt{4\pi*8.8\dots 10^{-12}\text{C}^2\text{N}^{-1}\text{m}^{-2}10^{-10}\text{m}10^3\text{J/mol}} \\ &=\sqrt{\frac{4/pi 8.8\dots 10^{-19}}{6.022\dots 10^{23}}}/1.602\cdot 10^{-19}=0.0268283 \vert e \vert \Rightarrow 1\vert e \vert = 37.274 \mathcal{Q},\end{align*}
    // $$1\nu=\mathcal{D}^3=\text{\AA}^3 \Rightarrow 1\text{\AA}^3 = 1\nu,$$
    // $$1\mathcal{C}=\mathcal{E}\mathcal{D}^6=1\text{kJ/mol}\text{\AA}^6 \Rightarrow 1\text{kJ/mol}\text{\AA}^6 = 1\mathcal{C},$$
    // where the permittivity of free space, $\epsilon_0$, was set to $8.854187817620\cdot 10^{-12}$C$^2$N$^{-1}$m$^{-2}$, the unit charge was set to $1.6021766208\cdot 10^{-19}$C and
    // we used a conversion factor 1Pa=$9.86923267*10^{-6}$atm. Note that all angles are in radians, both in real units and in the LJ like units.
    /////////////////////////////////////////////////////////////////////////////////////////////////
}

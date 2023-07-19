# Neutron Source Simulation
 The text is also available as .pdf: [neutronenbron.pdf](https://github.com/JBusink/Neutronenbron_Simulatie/files/12095379/neutronenbron.pdf).
## Introduction
In this document we describe a method for calculating the energy density and flux of neutrons in the water tank present at the 'natuurkundepracticum' at the Univeristeit van Amsterdam. A primary source of neutrons is located near the center of a water tank. This source contains a quantity of Am-241 (activity approximately 10 GBq) mixed with Be-9. The $\alpha$-particles emitted from the Am can be captured by the Be under formation of C-12 and a single neutron. The neutron flux from the source is $10^6$ per second. The neutrons are emitted with an energy in the range 3 to 10 MeV. The energy is subsequently moderated in collisions with the protons in the water molecules. Because the mass of the neutron is nearly identical to that of the proton, the energy is transferred very efficiently in collisions: a neutron colliding with a stationary proton loses, on average, half its kinetic energy in the collision. Neutrons at a larger distance from the center of the tank (the location of the neutron source) have a undergone more collisions. Consequently, we expect the energy distribution to be centered around lower energy as we move further away from the the source, towards the outside rim of the container. There are two reasons why knowledge of the energy distribution as a function of distance to the source is useful. The first is related to safety, we want the water tank to be large enough to ensure that the energy of the neutrons has been moderated down to thermal energies when they have diffused to the outer edge of the tank and subsequently escape to the environment. The second reason relates to the use of the setup in experiments. Targets can be placed in the tank at various distances from the source for activation experiments. The efficiency of this activation depends on the flux of impinging neutrons but also on (the distribution of) their energy. As both these quantities depend on position in the water tank it is important to calculate them. In the remainder of this document we describe the details of this calculation. The actual calculation takes the from of a Monte Carlo type simulation implemented as Julia code.


<details>
<summary> Description of elastic n-p collisions. </summary>

We shall make use of the fact that the masses of the neutron and proton are almost the same. Furthermore we realize that the protons are bound in the water molecules with energies of the order of just a few eV. This is very small compared to the MeV energies of the neutrons (at least initially). Hence we can approximate the protons as stationary targets. In addition we will use the fact that for all relevant energies the particles move at sub-relativistic speeds. Hence we may use the classical equation for the energy $E$ of each of the colliding particles:

$$E  = \frac{1}{2} m v^2, \tag{1}$$


with $m$ the mass of the neutron or proton and $v$ its corresponding velocity. 

The momentum $p$ of a particle can be related to its energy as:

$$p  = \sqrt{2\ E}, \tag{2}$$


where we have taken the mass of the particle as our mass unit for convenience.

To relate the momentum and energy of a neutron after a collision to the values before the collision we use the following method:

First we transform to the center of mass (cm) frame of the neutron and the target proton. In this frame the magnitude of the momentum of each of the colliding particles is $p/2$, half of the value of the impinging neutron. These momenta are oriented in opposite directions. This is illustrated in Fig. 1. An elastic collision corresponds to a random reorientation of the relative momentum of the two colliding particles in the cm-frame. To return to the laboratory frame we add $p/2$ to the component of the rotated momentum vector in the original direction of motion.  We are not interested in keeping track of the proton momentum after the collision. Hence, we discard it. The magnitude of the momentum of the neutron after the collision can take on any value between zero and its value before the collision. Hence, on average, the neutron loses half its energy. This implies that for a starting energy of several MeV, it takes roughly 20 collisions to reduce the neutron energy to thermal levels.

The above procedure is only valid as long as the energy the neutron is significantly larger than the binding energy of the proton in a water molecule. In this case the assumption that the proton is an isolated stationary target is justified. Ultimately the neutron will have lost enough energy to thermalize and from then on its energy can be obtained by randomly drawing from a thermal (Boltzmann) distribution. 
</details>


<details>
<summary> Steps in the Simulation</summary>

We now describe in a bit more detail the steps of the simulation. We want to calculate a simulated trajectory of neutrons which start at the point $(0,0,0)$ in our water tank with an energy of 5 MeV. We take the initial momentum to be ${\mathbf{p}}_0$ the vector $(p_x, p_y, p_z)$ with magnitude $p_0$ given by Eq. 2 and random direction. We then transform this momentum vector to the center of mass (cm) frame giving:


$$\mathbf{p}_{0 \, cm}= \frac{1}{2}(p_x, p_y, p_z). \tag{3}$$



A collision with a proton which is stationary in the lab-frame is described in the cm-frame as a rotation of the relative momentum vector of the neutron and the proton to a random direction. Consequently, the momentum the neutron after the collision in the cm-frame is given by:


<figure>
<img src= "https://raw.githubusercontent.com/JBusink/Neutronenbron_Simulatie/main/Figures/Figures_markdown/neutroncollision-1.png" alt="Trulli" style="width:100%">
<figcaption align = "center"><b>Fig.1 - Schematic representation of an elastic collision of a neutron an a stationary proton. The top left panel represent a neutron moving to the right with a momentum of magnitude p represented by the arrow. The panel in the top right is the same situation but viewed in the cm-frame. The bottom left panel represents the situation directly after the collision, again in the cm-frame. The magnitude of the momenta of the neutron and the proton in the cm-frame are unchanged in the collision; only the direction is changed. Finally we transform back to the lab-frame by adding p/2 to the horizontal components of the momentum of the neutron (blue) and the proton (red). We ignore the proton momentum. The new momentum of the neutron is used as initial momentum for the next collision. </b></figcaption>
</figure>

$$\mathbf{p}_{1,cm}= \frac{p_0}{2}(\sin \theta \cos \phi, \sin \theta \sin \phi, \cos \theta), \tag{4}$$

where the polar coordinates are drawn randomly from $\theta \, \epsilon \,\{0,\pi\}$ and  $\phi \,  \epsilon \, \{0,2 \pi\}$, respectively. 
Finally the new neutron momentum vector is taken back to the lab frame by adding ${\mathbf{p}}_0/2$ and renaming the resultant vector as ${\mathbf{p}}_1$. The proces can now be repeated in identical fashion starting with momentum ${\mathbf{p}}_1$ which results in a final momentum ${\mathbf{p}}_2$ and so on.

In addition to obtaining a series of momentum vectors we also want to keep track of the position where our simulated neutron ends up after each collision. We start out at the position of the neutron source which we take as the center ${\mathbf{r}}_0 = (0, 0, 0) $ of our coordinate system. The position ${\mathbf{r}}_1$ where the first collision takes place is calculated as follows:

$${\mathbf{r}}_1 = {\mathbf{r}}_0 +\hat{ {\mathbf{p}}}_0 \, d, \tag{5}$$


where $\hat{ {\mathbf{p}}}_0 = \frac{{\mathbf{p}}_0}{p_0}$ is a unit vector in the direction of the momentum  ${\mathbf{p}}_0$, and 
$d$ is the distance traveled to the position of the next collision which is drawn from the following probability distribution:

$$P(d) = \frac{1}{\lambda_{\mathrm{mf}} } \exp( - \frac{d}{\lambda_{\mathrm{mf}} }). \tag{6}$$

Here $\lambda_{\mathrm{mf}}$ is the mean free path which is related to the elastic collision cross section $\sigma$ and the number $n$ of target nuclei (i.e. protons) per unit volume through:

$$\lambda_{\mathrm{mf}} = \frac{1}{n\, \sigma}. \tag{7}$$


Again, as with the momentum calculation, the position after next collision can be calculated in an analogous way by replacing the indices $0$ and $1$ in Eq. 5 by $1$ and $2$, respectively.

The cross section $\sigma$ depends on the energy (and hence the momentum) of the neutron. We use the empirical approximate form:

$$\sigma (E) = \sigma_0 \left(\frac{E_0^2}{E_0^2 + E^2} +\frac{E_T^2}{ E^2} \right)^{\frac{1}{4}}, \tag{8}$$

where $E_0=4.5 \cdot 10^4$ eV is an upper threshold energy above which the cross section decreases, $E_T=0.038$ eV is a lower threshold at which the energy distribution is approximately thermal. Below $E_T$ the cross section increases with decreasing energy. For $E_T < E < E_0$ the cross section is approximately constant at a value of $\sigma_0 = 20$ barn. In fig. 2 a plot of the cross section versus energy is shown.

In simulating the trajectory of a neutron we keep track of its phase space coordinates (position and momentum) for a maximum of 500 collisions. The simulation is terminated earlier if the neutron  reaches a position which has a distance to the origin of 30 cm.  We take this to be the outer edge of the water tank. At this point we assume that the neutron escapes from the tank.The assumption of a spherical tank does not correspond to its actual shape but is simplifies  the simulation and we do no expect this simplification to greatly alter the outcome of the simulation. A second way in which we can lose the neutron is when it undergoes an inelastic collision according to the following reaction:

$$^1_0n + ^1_1H \rightarrow  ^2_1H + \gamma. \tag{9}$$

The cross section $\sigma_i$ for this neutron capture process is approximately given by:

$$\sigma_i (E) = \sigma_1\sqrt{\frac{E_1}{E}} ,\tag{9}$$

where $E_1=1$eV is a reference energy and $\sigma_1=4.3\cdot 10^{-2}$ barn is the cross section at that reference energy.
This inelastic process describes the dominant neutron removal process in water. The inelastic cross section is well below the elastic cross section for all energies but it is practically negligible at energies above about $10^4$ eV (six orders of magnitude smaller than the elastic cross section). At thermal energies and below the ratio of $\sigma_i$ and $\sigma$ is $1/83$, still small but not entirely negligible.

Ultimately we are interested in three quantities, all determined as a function of distance to the origin (we assume spherical symmetry throughout): the density, the energy density, and the flux per unit area. In particular the latter quantity is relevant for neutron activation experiments. The cross sections governing the inelastic processes depend on the neutron flux impinging on the target, as well as the energy of the neutrons.

In practical terms the density $\rho(r)$ at radius $r$ is determined as follows. First we divide the spherical volume which defines the interior of the water tank into $N$ (for example 30) discrete shells spanning from radius $r$ to  $r+\Delta r$. For each shell we have a counter to which we add 1 each time the neutron under scrutiny undergoes a collision inside this shell. Note that this can happen multiple times for each neutron as it is subject to a quasi-random walk which can return to the same shell more than once. This process is repeated cumulatively for many (for example one million) neutrons to build up sufficient statistics. The obtained number of observed collisions within each shell after $N$ particles, divided by $4 \pi r^2$, is proportional to the density. This `density' is clearly not normalized properly because it is proportional to $N$ and depends on $\Delta r$ for example, but for many applications this is not necessary. 

The second quantity we calculate is the flux $\phi$ through a surface. The numerical procedure for for this is the following: we use the same spherical shells as defined above in the calculation of the density $\rho(r)$. Every time a neutron passes through radius $r$ between to successive collisions a counter associated with radius r is raised by one. We define two such counters: one counts the number of particles which cross the shell from the inside to the outside (this is $\phi_i$) the other ($\phi_o$) counts the number of particles which cross the shell in the opposite direction, from outside to inside. We can define the net flux through the shell in two ways: $\phi_+ =\phi_i + \phi_o$ is the total number of particles which hit the shell, from either side. This quantity is relevant for activation experiment as we don't care from which side the neutron impinges on our target. The second quantity which we can define is $\phi_- =\phi_i - \phi_0$. This is the {\em net} flux through the shell. This quantity would be exactly equal to $N$ in the absence of absorption processes: each particle which is launched at the origin of the tank eventually passes once through every shell, unless is it removed by an absorption process. As it stands there is absorption through the inelastic process and, hence, $\phi_-$ is less the $N$.

The quantity $\phi_+$ can be converted into a flux density (flux per unit area) by dividing through $4 \pi r^2$. If we combine the known energy of the neutron each time it is located in the shell at $r$ with the obtained value for $\phi_+$ we obtain the flux as a function of energy. This is the quantity which we need to estimate activation efficiencies. 

<figure>
<img src= "https://github.com/JBusink/Neutronenbron_Simulatie/blob/main/Figures/Figures_markdown/Fig_cross_sections-1-1.png?raw=true" alt="Trulli" style="width:100%">
<figcaption align = "center"><b>Fig.2 - Left panel: cross section for elastic neutron-proton scattering in barn (blue), plotted versus neutron energy assuming the protons to be stationary targets (see Eq. 8). The yellow curve is the inelastic cross 
section for the reaction of Eq. 9 given by Eq. 10. The right panel shows the corresponding mean free paths in mm.</b></figcaption>
</figure>
</details>

<details>
<summary>Results
</summary>
</details>


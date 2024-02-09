Hummer established tested its method for ubiquitin(76aa). Our molecules are bigger: VEGF(134aa), TGF beta (112aa), PDGF (109aa)


VEGF, TGF beta, PDGF. Are they monomers?

Self-diffusion is a measure of the inherent mobility of a species in solution as a result of random Brownian motion and is quantified by the self-diffusion coefficient, D. Its value is influenced by molecular properties, particularly relating to interactions with the surrounding environment. For example, the self-diffusion coefficient of ethanol is lower than suggested by its molecular mass, due to the “stickiness” of the hydrogen bond associations. As such, the self-diffusion coefficient can provide insight into molecular interactions as well as the aggregation behaviour of pure components and their mixtures. Self-diffusion coefficients can be calculated directly through Green–Kubo relationships or by tracking the mean-squared displacement of a particle over a certain time interval166 and can be directly measured in experiments.167 An important source of error in these calculations is the neglect in considering the unexpected effect of the simulation box size168 and the effect of periodic boundary conditions.169 A distinction must be made between this property and the transport diffusivity, the proportionality constant relating a mass flux to the gradient (pressure, chemical potential, concentration, etc.) that induces it. These latter transport quantities require more refined calculation methods and/or the use of non-equilibrium simulations.170 [https://pubs.rsc.org/en/content/articlehtml/2023/cp/d2cp05423j]

"Tutorial" about gromacs (no Hummer unwrapping)[https://www.sciencedirect.com/science/article/pii/S0167732221018304]

Trajectorydata are typicallystoredin a wrappedform,withcoordinatesof the moleculesfoldedbackinto the centralMDbox.Hence,for the computationof the cms-MSDsof themolecules,the trajectorydata needto be properlyunwrapped.It has recentlybeenreportedthat certainunwrappingschemes,whenappliedto trajectoriesfromNpTsimulations,mightintroduceartifacts,whichleadto unphysicallylargeself-diffusioncoefficients.47To unfoldourNVTandNpTsimulationdata,we havereliedon a proceduresimilarto thetoroidal-view-preservingschemediscussedin ref 48(Bullerjahn,J. T.; vonBülow,S.; Heidari,M.; Hénin,J.;Hummer,G. UnwrappingNPTSimulationsto CalculateDiffusionCoefficients.J. Chem.TheoryComput.2023,19,   3406−3417.),   which prevents these effects.It is builtinto our homegrownsoftwarepackagesMOSCITOand MDorado.It is our experiencethatthe popularMDAnalysispackagedoesalso seemto handlethisissueproperly (Michaud-Agrawal,N.; Denning,E. J.; Woolf,T. B.; Beckstein,O. MDAnalysis:A toolkitfor the analysisof moleculardynamicssimulations.J. Comput.Chem.2011,32,   2319−2327). [https://pubs.acs.org/doi/full/10.1021/acs.jpcb.3c07540]

Website para descargar dinámicas de proteínas: 
[https://mmb.irbbarcelona.org/www/node/354]
[https://mmb.irbbarcelona.org/MoDEL/]
Ejemplo: buscar 1EPG y otros growth factors

Hummer tiene un estimador de la difusión pero no lo contrastó con la self-difusión de proterínas. Porqué? Hay que buscar datos reales de difusión.

Protein self diffusion? https://link.springer.com/article/10.1007/s12551-023-01108-y

Qui fanno la self diffusion di ~500 molecole (pure, solvente) con MD all atom https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.26975
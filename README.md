# LunarLargeCraterParameterDatabase
Lunar crater parameter database for large craters with rim diameters of 100 to 650 km, as well as python codes for gravity forward modeling and inversion calling pySHTOOLS. 
The preprint of the accompanying paper "Investigating the Influences of Crustal Thickness and Temperature on the Uplift of Mantle Material Beneath Large Impact Craters on the Moon" has been post on ESSOAr: https://www.essoar.org/doi/10.1002/essoar.10503183.1

Content of This Repository: 

CraterParameterStatistics7.xlsx (Main sheet) includes following crater parameters:
  - Name, Longitude, Latitude, Diameter (Losiak et al., 2009; Head et al., 2010) 
  - Location [FHT/SPA/PKT] (Jolliff et al., 2000)
  - CBA, CBA sd, Pvalue: Crater central Bouguer gravity measurements (Reference Model a) referenced to 60 km above local Moho. Mean crustal thickness is assumed to be 35 km. 
  - Inverted central Moho uplift (Reference Model a): inversion of the Moho relief from gravity was conducted using pySHTOOLS.
  - Mare existence, Mare percentage, Maxiumum possible thickness, Maxium possible CBA due to mare deposits (Nelson et al., 2014; Whitten and Head, 2015)
  - Crustal thickness, porosity, and grain density (Wieczorek et al., 2013) 
  - Crater age (Losiak et al., 2009)
   
Gravity Test Models: 
  - Model b: Use uniform crustal density for Bouguer correction
  - Model c: Apply crater diameter-dependent high-pass filter to crater Bouguer gravity
  - Model d: Reference to 1,738 km
  - Model e: Reference to 1 km above local topography
  - Model f: Reference to 70 km above local Moho and assume mean crustal thickness of 45 km 
  
Moho Uplift Test Models: 
  - Model b: Use uniform crustal density for inversion
  - Model f: Mean crustal thickness = 45 km
  - Model g: enhanced high-frequency filtering in the inversion 
  
CbaSHTools folder includes python codes to call pySHTOOLS (M. A. Wieczorek, M. Meschede, E. Sales de Andrade, I. Oshchepkov, B. Xu, and A. Walker (2018). SHTOOLS: Version 4.4, Zenodo, doi:10.5281/zenodo.592762) for gravity modeling and inversion: 
  - ForwardMare_VariableDensity.py: forward model the gravity attaction due to maximum possible mare deposits
  - FAA2Moho_UniformDensity.py: Invert for Moho relief from free-air gravity assuming uniform crustal density
  - FAA2Moho_LateralVariableDensity.py: Invert for Moho relief from free-air gravity assuming laterally variable crustal denisty model (Reference: https://github.com/MarkWieczorek/pyCrust)
  
  

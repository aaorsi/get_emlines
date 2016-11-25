# get_emlines

This code computes line luminosities *log(L [erg s-1])* based on an input *SFR [M_sun/yr]* and metallicity *Z*. Details on the physics and the implementation behind this can be found in [Orsi et al. 2014](http://adsabs.harvard.edu/abs/2014MNRAS.443..799O)

**Example of usage:**
```
from get_emlines import *
line = 'Halpha'
log_lum = get_emlines(line,sfr,Z)
```
You can access the list of lines available in, e.g.: 
`HIImodels/Lines/LineInfo_Levesque10`
By default, the list of lines available (and their names) are:

```
# Emission Lines contained in file LineData_Levesque10
# id	Emission Line			Central Wavelength
0  Lyalpha			1216.0
1  Hbeta			4861.0
2  Halpha			6563.0
3  OII_3727			3727.0
4  OII_3729			3729.0
5  OIII_5007			5007.0
6  OIII_4959			4959.0
7  OI_6300			6300.0
8  NII_6548			6548.0
9  NII_6584			6584.0
10  SII_6717			6717.0
11  SII_6731			6731.0
12  NeIII_3870			3870.0
13  CII_158um			1.5800
14 NII_205um 2.0500
```
The central wavelgnth of the last two FIR lines is in um, the rest in Angstroms.

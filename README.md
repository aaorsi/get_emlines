# get_emlines

This code computes line luminosities *log(L [erg s-1])* based on an input *SFR [M_sun/yr]* and metallicity *Z*. Details on the physics and the implementation behind this can be found in [Orsi et al. 2014](http://adsabs.harvard.edu/abs/2014MNRAS.443..799O)

**Quick example of usage:**
```python
import get_emlines as lines
from get_emlines import *
sfr  = 0.1   # Msun/yr
Z    = 0.01  # M_metals/M_gas 
lums = lines.get_emlines('x',sfr,Z,all_lines=True)

# OIII 5007 luminosity:
print lums['OIII_5007']

# All luminosities available should retrieve something like this:
lums.dtype
Out[3]: dtype([('Lyalpha', '<f4'), ('Hbeta', '<f4'), ('Halpha', '<f4'), ('OII_3727', '<f4'), ('OII_3729', '<f4'), 
('OIII_5007', '<f4'), ('OIII_4959', '<f4'), ('OI_6300', '<f4'), ('NII_6548', '<f4'), ('NII_6584', '<f4'), 
('SII_6717', '<f4'), ('SII_6731', '<f4'), ('NeIII_3870', '<f4'), ('CII_158um', '<f4'), ('NII_205um', '<f4')])

```
Note: In its current version, the first argument is obsolete but must be included (hence, `'x'`), and `all_lines` must be set to `True`. This is work in progress, so things will be cleaned up pretty soon.


**A quick installation:**

- First, clone the repository.
- Add it to your python repositories (e.g. in your bashrc:

`export PYTHONPATH="${PYTHONPATH}:/home/your-python-repositories/get_emlines/"`

- import it as in the example above.




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

**Modifying gamma and q0 from Orsi+14**

The ionization parameter is assumed to be related to the gas-phase metallicity by a power-law. You can change the slope `g0` or normalization `q0` from their standard values (from Orsi+14) by doing:
```python
loii_disk = lines.get_emlines('xxx',sfrdisk, zdisk, g0 = -1.5,q0 = 3.5e7,all_lines=True)
```
In the above the slope and normalization changed to `-1.5` and `3.5e7`, respectively.



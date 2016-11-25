# get_emlines

This code computes line luminosities *log(L [erg s-1])* based on an input *SFR [M_sun/yr]* and metallicity *Z*.

**Example usage:**
```
from get_emlines import *
line = 'Halpha'
log_lum = get_emlines(line,sfr,Z)
```

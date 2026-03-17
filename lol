3,4c3,4
< import globe
< from utils import (
---
> from . import globe
> from .utils import (
9,12c9,12
< from read_files import Read_LC_Coeffs, Read_Ugg, Read_Gas_Props
< from occupancy import Occupancies
< from phase_equilibrium import DissociationTemperature, DissociationPressure
< from output import (
---
> from .read_files import Read_LC_Coeffs, Read_Ugg, Read_Gas_Props
> from .occupancy import Occupancies
> from .phase_equilibrium import DissociationTemperature, DissociationPressure
> from .output import (
84c84
<     """Execute the chosen calculation, print results, and save to a CSV file."""
---
>     """Execute the chosen calculation, print results, and save to a file."""
211c211
<   python vdWP_FL.py
---
>   vdwpfl
214c214
<   python vdWP_FL.py --gas Methane --property 3 --T-range 273,290,10
---
>   vdwpfl --gas Methane --property 3 --T-range 273,290,10
217c217
<   python vdWP_FL.py --gas CO2 --property 0 --T 280 --P 5000000 --output my_results
---
>   vdwpfl --gas CO2 --property 0 --T 280 --P 5000000 --output my_results
251c251
<              "The .csv extension is added automatically if omitted.",
---
>              "The .dat extension is added automatically if omitted.",
258c258
< if __name__ == "__main__":
---
> def main():
268a269,272
> 
> 
> if __name__ == "__main__":
>     main()

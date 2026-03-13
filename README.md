# CAP2+

I present a Matlab-based enviroment, CAP2+ that uses concepts from order statistics to infer if a maximum magnitude is influencing a catalogue of earthquakes.  This repository is related to a study about the same topic [Schultz, 2026], building off of and improving an older idea [Schultz, 2024].  This repository contains programs, scripts, and data to recreate figures and results of the study.  

The core tools developed are KS_dM_test.m, M2fit.m, and EnsembleW.m.  See also my prior GitHub repositories for required routines (e.g., Bval.m).  Note that the NLEtable.mat and NLEtable_small.mat files will need to be regenerated - they are too large to share via GitHub, see script_MakeNLE_LookupTables.m to do so.

References: 

            R. Schultz, (2026)
            Improving the resolvability of Mmax truncation via deeper order statistics
            Geophysical Journal International, XX, xx.
            doi: 10.1093/gji/ggag110.
            
            R. Schultz, (2024)
            Inferring maximum magnitudes from the ordered sequence of large earthquakes
            Philosophical Transactions of the Royal Society A, 382, 20230185.
            doi: 10.1098/rsta.2023.0185.
            

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details: http://www.gnu.org/licenses/

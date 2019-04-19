
/*
    Solvent properties and utilities.

    Copyright (C) 2019 Simone Conti
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <thermo.h>

/*
    Solvents data.
    Taken from: A.J. Garza, "Solvation Entropy Made Simple", JCTC 2019
*/
ThermoSolvent thermo_solvents[] = {
{ .name="ethylene"       , .density=0.570   , .acentricity=0.089   , .permittivity=1.000   , .expansion=2.400   , .mass=28.053  , .bbox=101.392 , .vvdw=39.674  , .rgyr=1.315    },
{ .name="iodine"         , .density=3.960   , .acentricity=0.229   , .permittivity=4.000   , .expansion=0.000   , .mass=253.809 , .bbox=134.764 , .vvdw=58.904  , .rgyr=1.226    },
{ .name="cyclohexane"    , .density=0.778   , .acentricity=0.212   , .permittivity=2.000   , .expansion=1.210   , .mass=84.160  , .bbox=228.762 , .vvdw=100.710 , .rgyr=2.039    },
{ .name="benzene"        , .density=0.880   , .acentricity=0.212   , .permittivity=2.300   , .expansion=1.250   , .mass=78.112  , .bbox=199.623 , .vvdw=83.444  , .rgyr=2.013    },
{ .name="toluene"        , .density=0.867   , .acentricity=0.263   , .permittivity=2.400   , .expansion=1.080   , .mass=92.139  , .bbox=233.693 , .vvdw=100.262 , .rgyr=2.303    },
{ .name="m-xylene"       , .density=0.860   , .acentricity=0.325   , .permittivity=2.400   , .expansion=0.990   , .mass=106.166 , .bbox=264.478 , .vvdw=116.810 , .rgyr=2.568    },
{ .name="o-xylene"       , .density=0.880   , .acentricity=0.310   , .permittivity=2.600   , .expansion=0.000   , .mass=106.166 , .bbox=246.536 , .vvdw=116.701 , .rgyr=2.451    },
{ .name="pentane"        , .density=0.626   , .acentricity=0.251   , .permittivity=1.400   , .expansion=1.580   , .mass=72.149  , .bbox=208.747 , .vvdw=95.308  , .rgyr=2.325    },
{ .name="isopentane"     , .density=0.616   , .acentricity=0.227   , .permittivity=1.800   , .expansion=0.000   , .mass=72.149  , .bbox=229.972 , .vvdw=95.292  , .rgyr=2.100    },
{ .name="hexane"         , .density=0.660   , .acentricity=0.299   , .permittivity=1.900   , .expansion=1.410   , .mass=86.176  , .bbox=233.592 , .vvdw=112.064 , .rgyr=2.662    },
{ .name="octane"         , .density=0.703   , .acentricity=0.398   , .permittivity=2.000   , .expansion=1.140   , .mass=114.229 , .bbox=287.002 , .vvdw=145.610 , .rgyr=3.346    },
{ .name="chloroform"     , .density=1.490   , .acentricity=0.218   , .permittivity=4.800   , .expansion=1.270   , .mass=119.378 , .bbox=196.352 , .vvdw=75.822  , .rgyr=1.452    },
{ .name="dioxane"        , .density=0.796   , .acentricity=0.307   , .permittivity=2.300   , .expansion=1.120   , .mass=88.106  , .bbox=179.114 , .vvdw=85.513  , .rgyr=1.898    },
{ .name="acetaldehyde"   , .density=0.788   , .acentricity=0.303   , .permittivity=21.100  , .expansion=1.690   , .mass=44.053  , .bbox=126.436 , .vvdw=47.858  , .rgyr=1.390    },
{ .name="acetone"        , .density=0.784   , .acentricity=0.304   , .permittivity=20.700  , .expansion=1.430   , .mass=58.079  , .bbox=159.069 , .vvdw=64.454  , .rgyr=1.721    },
{ .name="ethyl_acetate"  , .density=0.810   , .acentricity=0.329   , .permittivity=6.000   , .expansion=1.380   , .mass=88.106  , .bbox=220.449 , .vvdw=89.664  , .rgyr=2.355    },
{ .name="acetic_acid"    , .density=1.050   , .acentricity=0.447   , .permittivity=6.200   , .expansion=1.100   , .mass=60.052  , .bbox=145.735 , .vvdw=55.674  , .rgyr=1.578    },
{ .name="acetonitrile"   , .density=0.786   , .acentricity=0.278   , .permittivity=37.500  , .expansion=1.360   , .mass=41.052  , .bbox=126.689 , .vvdw=45.678  , .rgyr=1.365    },
{ .name="dimethyl_ether" , .density=0.740   , .acentricity=0.200   , .permittivity=5.300   , .expansion=0.000   , .mass=46.069  , .bbox=145.924 , .vvdw=53.961  , .rgyr=1.609    },
{ .name="diethyl_ether"  , .density=0.713   , .acentricity=0.281   , .permittivity=4.300   , .expansion=1.600   , .mass=74.122  , .bbox=205.264 , .vvdw=87.340  , .rgyr=2.329    },
{ .name="helium"         , .density=0.130   , .acentricity=-0.365  , .permittivity=1.100   , .expansion=-1.490  , .mass=4.003   , .bbox=47.040  , .vvdw=11.511  , .rgyr=0.000    },
{ .name="neon"           , .density=1.210   , .acentricity=-0.029  , .permittivity=1.500   , .expansion=15.400  , .mass=20.180  , .bbox=56.918  , .vvdw=15.291  , .rgyr=0.000    },
{ .name="argon"          , .density=1.400   , .acentricity=0.001   , .permittivity=1.500   , .expansion=4.800   , .mass=39.948  , .bbox=84.826  , .vvdw=27.830  , .rgyr=0.000    },
{ .name="krypton"        , .density=2.410   , .acentricity=0.005   , .permittivity=1.700   , .expansion=0.000   , .mass=83.798  , .bbox=97.930  , .vvdw=34.464  , .rgyr=0.000    },
{ .name="xenon"          , .density=2.940   , .acentricity=0.008   , .permittivity=1.900   , .expansion=0.000   , .mass=131.294 , .bbox=111.974 , .vvdw=42.272  , .rgyr=0.000    },
{ .name="water"          , .density=1.000   , .acentricity=0.344   , .permittivity=78.500  , .expansion=0.210   , .mass=18.015  , .bbox=68.813  , .vvdw=19.413  , .rgyr=0.684    },
{ .name="methanol"       , .density=0.796   , .acentricity=0.556   , .permittivity=32.600  , .expansion=1.090   , .mass=32.042  , .bbox=107.024 , .vvdw=36.757  , .rgyr=1.192    },
{ .name="ethanol"        , .density=0.796   , .acentricity=0.644   , .permittivity=24.600  , .expansion=1.090   , .mass=46.069  , .bbox=144.540 , .vvdw=53.478  , .rgyr=1.558    },
{ .name="propanol"       , .density=0.803   , .acentricity=0.623   , .permittivity=20.100  , .expansion=0.790   , .mass=60.095  , .bbox=174.781 , .vvdw=70.182  , .rgyr=1.858    },
{ .name="isopropanol"    , .density=0.786   , .acentricity=0.665   , .permittivity=17.900  , .expansion=0.000   , .mass=60.095  , .bbox=174.969 , .vvdw=70.064  , .rgyr=1.763    },
{ .name="butanol"        , .density=0.810   , .acentricity=0.593   , .permittivity=17.800  , .expansion=0.750   , .mass=74.122  , .bbox=199.101 , .vvdw=86.981  , .rgyr=2.179    },
{ .name="isobutanol"     , .density=0.802   , .acentricity=0.592   , .permittivity=17.300  , .expansion=0.940   , .mass=74.122  , .bbox=214.029 , .vvdw=87.020  , .rgyr=1.992    },
};

/*
    Return the id of a solvent from the solvent name.
    Returns -1 if failed.
*/
int thermo_get_solvent_from_name(char *name) {
    int nsolvents = sizeof(thermo_solvents)/sizeof(ThermoSolvent);
    int i;
    for (i=0; i<nsolvents; i++) {
        if (strcmp(thermo_solvents[i].name, name)==0) {
            return i;
        }
    }
    return -1;
}

/*
    Get solvent data from id.
*/
ThermoSolvent thermo_get_solvent_from_id(int id) {
    return thermo_solvents[id];
}

/*
    Print to fp information about the solvent id.
*/
void thermo_print_solvent_id(FILE *fp, int id) {
    if (id<0 || id>(int)(sizeof(thermo_solvents)/sizeof(ThermoSolvent))) {
        fprintf(fp, "ERROR! No solvent corresponds to ID %d\n", id);
        return;
    }
    ThermoSolvent solvent = thermo_solvents[id];
    thermo_print_solvent(fp, solvent);
    return;
}

/*
    Print to fp information about the given solvent.
*/
void thermo_print_solvent(FILE *fp, ThermoSolvent solvent) {
    fprintf(fp, "Solvent:                 %s\n", solvent.name);
    fprintf(fp, "   Molecular weight:     %g g/mol\n", solvent.mass);
    fprintf(fp, "   Density:              %g g/ml\n", solvent.density);
    fprintf(fp, "   Acentric factor:      %g\n", solvent.acentricity);
    fprintf(fp, "   Permittivity:         %g\n", solvent.permittivity);
    fprintf(fp, "   Isobaric expansion:   %g *1E-3 /K\n", solvent.expansion);
    fprintf(fp, "   Gyration radius:      %g A\n", solvent.rgyr);
    fprintf(fp, "   Bounding box area:    %g A^2\n", solvent.bbox);
    fprintf(fp, "   Van der Waals volume: %g A^3\n", solvent.vvdw);
    fprintf(fp, "\n");
    return;
}


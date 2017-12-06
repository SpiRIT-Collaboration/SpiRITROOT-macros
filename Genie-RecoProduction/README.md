# create_submit.C

This macro creates submit files for NSCL cluster with the correct gating grid noise data file and pulser calibration file.

**You must modify some part of the macro to adapt to your system.**

You can search `Fix me` inside the macro to find where to modify.

`run_info.txt` file contains run numbers and corresponding total event numbers. It should be in the same folder as `create_submit.C` macro file.

# makeBxBz0.sh

This script reads in SAMURAI magnetic field map file and outputs the field map file having Bx=Bz=0.

# run_reco_experiment.C

This macro file is used to generate reconstruction ROOT file production.

To get better momentum resolution, the following requirement should be satisfied.
- All layers are used. (112 layers)
- Bx=Bz=0 field map file is used.
- Magnetic field map should be shifted to upstream by 25 cm in (STGFBField.cc)[https://github.com/SpiRIT-Collaboration/SpiRITROOT/blob/develop/reco/Genfit/STGFBField.cc]. (58 -> 33)

# RF_MLOC

Radiative forcing (RF) model applicable to multiple cloud layers overlapping (Multiple-Layer Overlapping Clouds)

This radiative forcing model has been developed with the objective of evaluating the impact on contrail-attributable radiative forcing of multiple layers overlap. This includes the evaluation of the effect of contrail-contrail overlaps on contrail RF and the effect of cloud-contrail overlaps on contrail RF.

The model has been verified with state-of-art clouds RF models.

This code was developed in 2019 by Inés Sanz-Morère, grad student from the Laboratory of Aviation and the Environment (LAE) at MIT (inessanz@mit.edu).

Input data required
--------------

- For each contrail and cloud layer we require the following information:
1. Vertical level (or relative location with respect to other layers) (from bottom to top)
2. Area (in m2)
3. Ambient temperature at that level, T (in K)
4. Optical depth, od
5. Asymmetry parameter, g (representing microphysical properties)

This information should be saved in a structure with the number of elements equal to the number of layers.

- As atmospheric data we require:
1. Outgoing longwave radiation, OLR (in W/m2)
2. Solar direct radiation, SDR (in W/m2)
3. Atmospheric transmittance of solar radiation, t
4. Surface albedo, alpha
5. Cosinus of solar zenith angle, mu


Output data information
--------------

Based on the simple cloud radiative forcing model developed by Corti and Peter (2009), with atmospheric and clouds data, and contrail characteristics, the model is capable of providing four different RF outputs related to contrail impact:
1. Only clouds RF (in W)
2. Only contrails (clear sky) (in W)
3. All layers (clouds and contrails) (in W)
4. Independent contrails (as if they were not overlapping) in all-sky (in W)

Output data are under the form of a structure of four elements with the following three fields: longwave RF (in W), shortwave RF (in W) and maximum area (in m2) of all layers considered.

Substracting output 1 to output 3 provides only contrails RF in an all sky situation (with the presence of natural clouds), accounting for them overlapping. 
The difference between this last value and output 4 informs on the impact of contrail-contrail overlap on contrail RF.
The difference between this last value and output 2 informs on the impact of cloud-contrail overlap on contrail RF.


Usage
--------------

Download the script and run it with the required input data.
If RF in W/m2 is looked for, all areas should be set to 1.

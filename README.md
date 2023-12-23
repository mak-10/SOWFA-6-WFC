# SOWFA-6-WFC
Includes a wind farm canopy as an fvOption to src/ABLforcing/windCanopy.
A WFC uniformly distributes the cumulative thrust force of a wind farm over its volume.
Requires toposet to choice the wind farm or canopy control volume as a cellSet. Thus, relies on the class cellSetOptions.
The constant/fcOptions file gives the Cft, which is the wind farm thrust coefficient, canopyArea, and direction of the thrust force. 

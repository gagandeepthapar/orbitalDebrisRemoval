# **AERO - 351 Final Project**
**Orbital Mechanics:** Hypothetical orbital debris rendezvous demonstration using propulsive maneuvers.

## **Background**
AERO-351, or Orbital Mechanics, is a class taught at the California Polytechnic State University at San Luis Obispo (Cal Poly SLO) which serves as an introduction to...
> "[The] Motion of a body in a central field. Keplerian Orbits. Orbital Manuevers. Interplanetary trajectories." - Cal Poly Aerospace Engineering Course Catalog

The final project in the course was a group effort to select 4 pieces of orbital debris from an online catalog, rendezvous with each debris piece, and determine the Delta-V required to complete the mission. No specific constraints were given as the final project was a demonstration of the material taught in the course, however, orbit progression for each of the debris must be accounted for before determining the rendezvous location and Delta-V required.

The repository contains the MATLAB files for the project and associated figures.

## **Criteria:**</br>
The final project required a set of criteria to be achieved to be considered successful:
* Select 4 pieces of debris; two must be in LEO, one in MEO, and one in GEO
* Debris objects *must* be inactive/dead spacecraft/debris; i.e., they cannot be active satellites
* The departure date can be arbitrary
* The GEO object must have less than 5 deg of inclination but greater than 0 deg
* Once the spacecraft has rendezvous'd with its target debris, it must stay with the debris for at least 5 periods
* It can be assumed there are no orbital perturbations for any of the debris
* The total Delta-V must be less than **18 km/s**

## **TLE Positional Data (2nd Line) for...**</br>
### **LEO Debris (#1): Iridium-33 Debris:**</br>
``[24946 | 86.3843 | 127.9418 | 0008700 | 151.1544 | 209.0134 | 14.33702974211628]``

### **LEO Debris (#2): Iridium-33 Debris:**</br>
``[33776 | 86.4036 | 138.4324 | 0015334 | 156.4007 | 214.2811 | 14.34129899613840]``

### **MEO Debris: GLONAS Rocket Body:**</br>
``[13610 | 64.0303 | 137.2256 | 0008118 | 199.9412 | 343.5184 | 2.14005188297957]``

### **GEO Debris: INTELSAT 2-F2:**</br>
``[02639 | 1.9465 | 287.8917 | 0009103 | 316.2065  | 67.8614 | 1.0031297298739]``


## **Orbital Manuevers (Start Time JD 21 181.78619214501):**</br>
### **LEO Debris (#1) to LEO Debris (#2):**</br>
* Circularize at LEO #1 apogee
* Perform impulse burn to move into LEO #2 orbital plane (inclination and RAAN change)
* Simultaenously initiate Hohmann Transfer to end at LEO #2 apogee
* Decircularize into LEO #2 when the debris reaches its apogee (accounts for change in Arugment of Perigee)
* Perform phasing maneuver to account for change in True Anomaly to rendezvous with debris
* Stay "attached" to debris in LEO #2 orbit for 5 periods (as required by the project)</br></br>
**DELTA-V: 3.1825 km/s**</br>
**DELTA-T: 15.0614 hrs**</br>
### **LEO Debris(#2) to MEO Debris </br>**

* Two-Impulse Burn to chase MEO Debris via Lambert's Solution
* Stay "attached" to debris in LEO #2 orbit for 5 periods (as required by the project)</br></br>
**DELTA-V: 4.4022 km/s**</br>
**DELTA-T: 4.7917 days**</br>
### **MEO Debris to GEO Debris </br>**

* Two-Impulse Burn to chase GEO Debris via Lambert's Solution
* Stay "attached" to debris in LEO #2 orbit for 5 periods (as required by the project)</br></br>
**DELTA-V: 9.1105 km/s**</br>
**DELTA-T: 21.76 days**</br>

## **TOTAL TIME and DELTA-V REQ'D:**
``TOTAL DELTA-V: 16.6952 km/s``</br>
``TOTAL DELTA-T: 26.4067 day``

## **Potential Optimizations:**
The current setup only accounts for the earliest time the spacecraft can perform a burn to move to its next target. That is, as soon as its 5 periods with each debris is complete, the spacecraft initiates the burn to move on. Instead, had the orbits been propagated further than 5 periods, a more optimal path could be discovered. This is similar to what current space missions do; NASA doesn't launch their Mars rovers as soon as its complete. They wait until the best Mars pass considering the time it will take for the launch vehicle/spacecraft to reach Mars will be.

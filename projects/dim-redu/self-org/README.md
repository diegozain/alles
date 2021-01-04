# Self-organizing maps
diego domenzain
September 2020 @ Colorado School of Mines

## Cluster n-dimensional data in two-dimensions

Given a data set in n-dimensions, how do we visualize it in our limited human vision?

__These scripts are examples of self-organizing maps.__

* ```synth_ex.m``` is a light-weight synthetic example.
 
* ```zodiac_ex.m``` uses data from one of those buzz-feed-type questionnaires but for a work environment.

  * Columns are the questions, rows are people.
  * Interesting trends:
    * Persons 1, 9, and 14 were the oldest ladies in the group.
    * Persons 2, 3, and 8 were officemates (all dudes).
    * Persons 5, 6, 7, and 11 were the most boring in the group by far (ladies and dudes).
    * Persons 10, 12, and 13 were international students (all dudes).
    
* ```geoclasses_ex.m``` uses data from _geoscience_ students (not necessarily geo-__physicists__). 
 
  * Columns are _classes_ they take, and rows are _attributes_ from those classes.
  * Interesting trends
    * Classes 2 and 12 are _hydrology_ and _geophysics_.
    * Classes 1, 3, 8, and 14 are _sed/strat_, _geomorphology_, _mineralogy_ and _earth materials_.
    * Classes 4 and 6 are _structure and field_.
    * Classes 9 and 10 are _historical_ and _paleontology_.
    * Classes 11 and 13 are _environment_ and _climate_.

---

## Synthetic example

[![](../pics/self-org-synth.png)](./)

## Zodiac example

[![](../pics/self-org-zodiac.png)](./)

## Geo-classes example

[![](../pics/self-org-geoclasses.png)](./)


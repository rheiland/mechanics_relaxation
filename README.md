# PhysiCell: an Open Source Physics-Based Cell Simulator for 3-D Multicellular Systems

```
(base) M1P~/git/mechanics_relaxation$ make load PROJ=cells11_final
(base) M1P~/git/mechanics_relaxation$ make 
(base) M1P~/git/mechanics_relaxation$ project  
or:
(base) M1P~/git/mechanics_relaxation$ pcstudio  

(base) M1P~/git/mechanics_relaxation$ python analysis/plot_11cells_crop.py 88
```
<img src="images/physicell_relax11.png" width=400/>

```
- assuming we have Chaste's results
(base) M1P~/git/mechanics_relaxation$ python analysis/plot_11cells_csv.py pc_plot_11cells.csv
```
<img src="images/physicell_vs_chaste.png" width=500/>

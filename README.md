# Using Inverse Optimization to Learn Cost Functions in Generalized Nash Games
## README

This code repository is broken down into files that require the 
following programs:
- MATLAB
- PYTHON
- GAMS

It also assumes that we are running on a Ubuntu/Linux operating system. 

The IO\_GNEP\_* shell scripts will run different collections of the other scripts
housed in this repository.  The bash scripts correspond to the four
experimental groups discussed in the paper. The Python classes called 
by the Python scripts are:

- RK\_IO\_methods which is housed in the RK\_IO\_model.py file.
This class contains the inverse optimization residual model for
the paper for both the same costs across all players and different
costs across players
- generalized\_RK\_framework which is housed in the Generalized\_RK\_Framework.py file
This class contains the mechanics methods for translating between
Python and the various other programs.  It also runs the 
RK\_IO\_methods, calculates total flow error, and saves all of the data.

**Note:** We have three class objects called class\_object\_1 through 3
that are created with the various scripts that call methods from
these classes.  They transfer data between the different Python scripts.

---

The Python scripts are named in terms of steps, so 
script\_1\_generating\_data\_for\_matlab.py indicates the first step.
However, there are MATLAB and GAMS calls that happen between
the various Python scripts.  For a general experiment, this is the
sequence of calls and their outputs in the bash scripts:

1. python script\_1\_generating\_data\_for\_matlab.py [with potentially arguments]
OR python script\_1\_generating\_data\_for\_matlab\_sioux\_falls.py
OUTPUT: Pickle object (class\_object\_1) and MATLAB data file
(data\_for\_matlab\_original\_costs.mat)

2. matlab script call to either matlab\_generating\_data\_original\_costs.m
or matlab\_generating\_data\_original\_costs\_different\_costs.m
OUTPUT: GDX files traffic\_data\_gdx\_iteration\_1.gdx through 10

3. gams script call 10 times for 10 trials, either to the 
gams\_original\_costs\_file.gms or gams\_original\_costs\_file\_different\_costs.gms
OUTPUT: CSV files traffic\_results\_original\_costs\_1.csv through 10

4. python script call to either script\_2\_inverse\_optimization\_step.py
or script\_2\_inverse\_optimization\_step\_different\_costs.py
OUTPUT: Pickle object (class\_object\_2)

5. matlab script call to either matlab\_generating\_data\_IO\_costs.m
or matlab\_generating\_data\_IO\_costs\_different\_costs.m
OUTPUT: GDX files IO\_traffic\_data\_gdx\_iteration1.gdx through 10

6. gams script call 10 times for 10 trials, either to the 
gams\_IO\_costs\_file.gms or gams\_IO\_costs\_file\_different\_costs.gms
OUTPUT: CSV files IO\_traffic\_results\_1.csv through 10

7. python script call to script\_3\_calculating\_flow\_error\_and\_boxplots.py
OUTPUT: Named CSV file of the data collected over the course
of the experiment and Pickle object (class\_object\_3) 

---

Note that the data from the experiments we ran can be found in the
csv\_data\_files\_from\_experiments folder.  The scripts in that folder
will also produce the boxplots for the 4 experimental groupings we 
describe in the paper.  We do use this opportunity to tell the reader
that the data for data\_saving\_experiment\_different\_costs\_5x5\_Grid\_num\_players\_10\_alpha\_5.csv
was inputted by hand due to problems with that experiment (as discussed
in the paper)

**ALSO:** The csv data files in this folder (along with the .mat and .gdx data
files) relate to the following run: 
Sioux Falls, number of players N=10, and alpha=10.  If a user would like to run experiments with the Sioux Falls network, they will need to download the `SiouxFall_flow.tntp` file separately from the [https://github.com/bstabler/TransportationNetworks](https://github.com/bstabler/TransportationNetworks "TransportationNetworks") repository.

See `1_PACKAGES.txt` for the Packages utilized by the two machines we utilized.

---

## Major References

Ratliff, Lillian J., et al. "Social game for building energy efficiency: Incentive design." 2014 *52nd Annual Allerton Conference on Communication, Control, and Computing (Allerton)*. IEEE, 2014.

Konstantakopoulos, Ioannis C., et al. "A robust utility learning framework via inverse optimization." *IEEE Transactions on Control Systems Technology* 26.3 (2017): 954-970.

Transportation Networks for Research Core Team. *Transportation Networks for Research*. https://github.com/bstabler/TransportationNetworks. Accessed 2020.
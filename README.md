# ne8abs_paper
code and scripts for finding main halos, making histograms and 2D
projected maps, and analysing those products in FIRE simulations.
This is a subset of the scripts in my fire_an repo which I used
specifically for Wijers et al. (2024), copied on July 23, 2024 from
the fire_an main branch.



'installation':
---------------
This is a selection of scripts, not a polished package. Many scripts, especially those in `queuerun` and 
`makeplots`, contain functions and other bits of code that I did not
end up using in the paper.

I have provided
a `condaenv_firep1.yml` file that includes the python packages required 
to run these scripts. It can be used to set up a conda environment 
where this code should run.

Note that a few functions requires additional packages, such as Andrew
Wetzel's FIRE analysis code. I haven't included these here, since these
are no longer the defaults I use.

Note that the main projection code and ion balance calculations call 
some compiled C code. I don't have that publicly available right now, 
because I did not write the SPH projection code myself (I only made 
small modifications), and I have not gotten around to asking the authors
for permission.

setup:
------
The `utils` directory should contain a file called `opts_locs.py`. This
file simply lists a number of files and directories needed to run the
code, such as where to search for simulation outputs, and where to store
halo data. Some example files `opts_locs_<system>.py` are provided.
Note that the `kernel_list`, `desngb`, `solar_abunds_ea`, and `Zsun_ea`
variables do not need to be changed. They're basically just stored here
because this started off as a parameter file.

Storing files as `opts_locs_<system>.py` and copying those to 
`opts_locs.py` locally is recommended, as it provides a back-up for the 
setups on github, without constantly overwriting local locations. 
(`utils/opts_locs.py` is included in .gitignore for this reason.)

running the scripts:
--------------------
To run, for example, `run.py`, options are:
- run from the directory above `ne8abs_paper` as
  ```
  python -m ne8abs_paper/queuerun/run -3
  ```
  (This will raise an error because `run.py` only takes postive integer
  arguments, but it also won't accidentally start an hour-long process
  on a login node.)
- or add the directory above ne8abs_paper (e.g. `foo`) to `PYTHONPATH`:
  ```
  export PYTHONPATH="/path/to/foo/:${PYTHONPATH}"
  ```
  and then, from anywhere, run: 
  ```
  python -m ne8abs_paper.queuerun.run -3
  ```
  or specify an absolute or relative path to `run.py` from your 
  current directory:
  ```
  python /path/to/ne8abs_paper/queuerun/run.py -3
  ```

This is basically an issue of the different scripts in different 
directories needing to be able to find each other. Running from the
directory above `ne8abs_paper` adds that directory to `PYTHONPATH`, like 
adding it to `PYTHONPATH` explicitly does. 

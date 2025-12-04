I report here the set of commands to run (Snakefile will need an update).

Create workspace:

```
bin/Lb2pktaul_dataset_build_workspace
```

Run Prob method:

```
bin/Lb2pktaul_dataset --var branchingRatio --npoints 100 --scanrange 0.:4.e-5 --CL 68.3 --CL 95.4 [--cls 2]
```

`[--cls 1]` and `[--cls 2]` give the same result with Prob method, but it is not needed for the plugin method (while Prob result is needed).

Run plugin on batch condor:

```
bin/Lb2pktaul_dataset -a pluginbatch  --var branchingRatio --npoints 100 --scanrange 0.:4.e-5 --ntoys 20 --nbatchjobs 100 --batchstartn 1 --batchreqs ../scripts/cern_condor_req.txt --batchsubmit
```

For Lb2pktaul (i.e. low BDT cut), 1 hour is enough so also `--batchreqs ../scripts/cern_condor_req_1hour.txt` can be used.

Equivalent running in local is:

```
bin/Lb2pktaul_dataset -a pluginbatch  --var branchingRatio --npoints 100 --scanrange 0.:4.e-5 --ntoys 20  --nrun 1
bin/Lb2pktaul_dataset -a pluginbatch  --var branchingRatio --npoints 100 --scanrange 0.:4.e-5 --ntoys 20  --nrun 2
...
bin/Lb2pktaul_dataset -a pluginbatch  --var branchingRatio --npoints 100 --scanrange 0.:4.e-5 --ntoys 20  --nrun 100
```

Read plugin results:

```
bin/Lb2pktaul_dataset -a plugin --var branchingRatio --npoints 100 --scanrange 0.:4.e-5 -j 1-100  --CL 68.3 --CL 95.4  --cls 2
```

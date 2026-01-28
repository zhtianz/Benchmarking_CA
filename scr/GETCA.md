# Methods for Obtaining $C\alpha$ Atom Prediction Results

This document describes how to extract $C\alpha$ atom prediction results from various atomic model building methods. Please note that for benchmarking purposes, the test files used are the results obtained after $C\alpha$ prediction and subsequent clustering.

## 1. ModelAngelo
**Repository:** [https://github.com/3dem/model-angelo](https://github.com/3dem/model-angelo)

To only run the $C\alpha$ atom prediction and skip the post-processing refinement, modify the configuration file before execution.

- **Modification:** 
  In `nucleotides_no_seq/config.json`, change `"gnn_infer_args": {"num_rounds": 3}` to `"num_rounds": 0`.
- **Execution:**
  ```bash
  conda activate model_angelo
  model_angelo build_no_seq -v [map_path] -o [output_path] --keep-intermediate-results
  ```
- **Output File:** `${output_path}/see_alpha_output/see_alpha_output_ca.cif`

---

## 2. CryoAtom
**Repository:** [https://github.com/SBQ-1999/CryoAtom](https://github.com/SBQ-1999/CryoAtom)

Similar to ModelAngelo, the iterative rounds must be disabled to extract the initial $C\alpha$ predictions.

- **Modification:** 
  In `CryoAtom/config.json`, change `"CryNet_args": {"num_rounds": 3}` to `"num_rounds": 0`.
- **Execution:**
  ```bash
  conda activate CryoAtom
  cryoatom build --map-path [map_path] --fasta-path [fasta_path] --output-dir [output_path] --keep-intermediate-results
  ```
- **Output File:** `${output_path}/see_alpha_output/see_alpha_output_ca.cif`

---

## 3. DeepMainMast
**Repository:** [https://github.com/kiharalab/DeepMainMast](https://github.com/kiharalab/DeepMainMast)

Truncate the process after the predicted $C\alpha$ clustering in the original code, then run the following command:

- **Execution:**
  ```bash
  conda activate deepmainmast
  ./dmm_full_multithreads.sh -p [program_path] -m [map_path] -f [fasta_path] -c [contour] -o [output_path] -t [path_training_time] -T [fragment_assembling_time] -C [num_cpu] -M [num_cpu]
  ```
- **Output File:** `${output_path}/results/NODE_p0.3.pdb`

---

## 4. EModelX
**Repository:** [https://bio-web1.nscc-gz.cn/app/EModelX](https://bio-web1.nscc-gz.cn/app/EModelX)

Truncate the process after the predicted $C\alpha$ clustering in the original code, then run the following command:

- **Execution:**
  ```bash
  conda activate EModelX
  python run.py --protocol=temp_free --EM_map=[map_path] --fasta=[fasta_path] --output_dir=[output_path] --run_pulchra --pulchra_path modules/pulchra304/src/pulchra
  ```
- **Extraction Method:** 
  The $C\alpha$ atom results are extracted from the point set generated after the `clustering()` function and saved into a new PDB/CIF file.

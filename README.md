# Benchmarking Deep Learning Methods for Cα Atom Prediction in Cryo-EM Density Maps

## Project Overview
This project aims to evaluate the performance of Cα atom prediction across different atomic model building methods.

---

## Installation and Environment Setup

### Dependencies
| Dependency | Version |
|------------|---------|
| Python     | 3.8.8   |
| Numpy      | 1.24.1  |
| Scipy      | 1.10.1  |
| torch      | 1.13.1  |
| Bio        | 1.78    |
| einops     | 0.8.0   |
| tqdm       | 4.66.5  |

### Installation Steps
   ```bash
   cd <your_project_path>
   git clone https://github.com/zhtianz/Benchmarking_CA.git
   ```

## Usage
   ```bash
   cd Benchmarking_CA/src
   python CAEVAL.py -r <native_model.pdb> -p <predict_model.pdb> -o <output_file>
   ```

### Commands
- `-r` : Reference model in PDB/CIF format
- `-p` : Predicted model in PDB/CIF format
- `-o` : Output file

## Example
   ```bash
   cd Benchmarking_CA/src
   python CAEVAL.py -r example/7kjr.pdb -p example/22898_pred.cif -o eval_result.log
   ```

The evaluation results can be found in the file  `eval_result.log`.

---

## Additional Information

### Obtaining Prediction Results
For detailed instructions and methods on how to obtain $C\alpha$ atom prediction results, please refer to the documentation in [src/GETCA.md](src/GETCA.md).

### Datasets
The `datasets` directory contains comprehensive details regarding the **test datasets** used in this benchmark.


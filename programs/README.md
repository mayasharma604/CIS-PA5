# CISPA-5
Computer Integrated Surgery Programming Assignment 5

**Authors:** Maya Sharma
**Date:** 12/6/25

## Installation Instructions
1. Make sure Python 3.8 or higher is installed on your system.
2. Install required dependencies:
   ```bash
   pip install numpy
   ```
3. Place the input data files (`Problem5-BodyA.txt`, `Problem5-BodyB.txt`, `Problem5MeshFile.sur`, `Problem5Modes.txt`, `PA5-X-{Debug/Unknown}-SampleReadingsTest.txt`) in the `2025_PA345_Student_Data` folder relative to the source files.

## File Descriptions
**pa5.py** is the main script that runs the deformable registration algorithm for a single dataset. It accepts a file letter (A-F for debug, G/H/J/K for unknown) and processes the corresponding sample readings file.

**deform_registration.py** contains the core deformable registration implementation:
- `compute_barycentric()`: Computes barycentric coordinates of a point within a triangle
- `read_modes_fixed()`: Parses statistical shape model files with mean shape and deformation modes
- `solve_pa5()`: Main deformable registration algorithm combining rigid and non-rigid transformations

**utility_functions.py** contains helper functions including:
- File parsing for bodies, mesh, and sample readings
- Point-to-point registration (Arun's method)
- Frame transformations (`apply_transform`, `apply_inverse_transform`)
- Mesh query functions (`closest_point_on_triangle`, `closest_point_on_mesh`)

**ICP_algo.py** contains PA4 ICP implementation reused for initial rigid registration and core algorithms.

**unit_tests.py** contains comprehensive tests validating triangle projection, mesh queries, and algorithm correctness.

## Execution Instructions
To run the program for a specific dataset:
```bash
python pa5.py <letter>
```
Examples:
```bash
python pa5.py A    
python pa5.py G    
```

**Valid letters:** A-F (debug datasets), G, H, J, K (unknown datasets)

The output file `PA5-X-{Debug/Unknown}-Output.txt` will be created containing:
- Header: `N_samps PA5-A-Debug-Answer.txt N_modes`
- Mode weights: `λ₁ λ₂ ... λ_Nmodes` (4 decimal places)
- For each sample: `s_x s_y s_z    c_x c_y c_z    |s_k - c_k|`

## Dependencies
- **NumPy** (>=1.20.0): For numerical operations and linear algebra
- **Standard Python libraries**: sys, os, re


## References
```bibtex
@misc{benjamindkilleen2022Sep,
 author = {Killeen, Benjamin D.},
 title = {{cispa: Template for CIS I programming assignments at Johns Hopkins}},
 journal = {GitHub},
 year = {2022},
 month = {Sep},
 url = {https://github.com/benjamindkilleen/cispa}
}
```

## Notes
- The program implements deformable registration combining rigid ICP with statistical shape model deformation
- Convergence tolerance defaults to 0.001 with maximum 50 iterations
- Output format matches PA5 requirements with 4 decimal places for λ values
- All mathematical derivations follow the assignment specifications

Killeen, B. D. (Sept. 8, 2022). Frame Transformations. Retrieved from https://benjamindkilleen.com/files/frame_transformations.pdf

K. S. Arun, T. S. Huang, and S. D. Blostein, "Least-squares fitting of two 3-D point sets," IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. PAMI-9, no. 5, pp. 698–700, Sept. 1987.

Taylor, R. (1-73). Medical Robots Part 2 [Lecture presentation]. Computer Integrated Surgery 1, Johns Hopkins University

Taylor, R. (2025). Notes for PA5: Deformable Registration to a Statistical Shape Model. Johns Hopkins University.

Cootes, T. F., & Taylor, C. J. (2000). Statistical Models of Appearance for Computer Vision.
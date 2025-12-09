# pa5_single_run.py
import os
import sys
from deform_registration import solve_pa5_correct


if len(sys.argv) != 2:
    print("Usage: python pa5.py <File_Letter>")
    print("Example: python pa5.py A")
    print("Example: python pa5.py G")
    sys.exit(1)

file_letter = sys.argv[1].upper()
if file_letter not in "ABCDEFGHJK":
    print(f"Error: Invalid file letter '{file_letter}'. Use A-F for debug or G, H, J, K for unknown.")
    sys.exit(1)

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir) 

data_folder = os.path.join(parent_dir, "2025_PA345_Student_Data")
if not os.path.exists(data_folder):
    data_folder = os.path.join(parent_dir, "2025 345 Student Data") 
    
if not os.path.exists(data_folder):
    print(f"ERROR: Could not find data folder.")
    sys.exit(1)

is_debug = file_letter in "ABCDEF"
file_type = "Debug" if is_debug else "Unknown"

sample_filename = f"PA5-{file_letter}-{file_type}-SampleReadingsTest.txt"
output_filename = f"PA5-{file_letter}-{file_type}-Output.txt"

sample_file_path = os.path.join(data_folder, sample_filename)
output_file_path = os.path.join(current_dir, output_filename)

bodyA_file = os.path.join(data_folder, "Problem5-BodyA.txt")
bodyB_file = os.path.join(data_folder, "Problem5-BodyB.txt")
mesh_file = os.path.join(data_folder, "Problem5MeshFile.sur")
modes_file = os.path.join(data_folder, "Problem5Modes.txt")

if not os.path.exists(sample_file_path):
    print(f"ERROR: Sample file not found: {sample_file_path}")
    sys.exit(1)


try:
    s_k, c_k, lambdas = solve_pa5_correct(
        bodyA_file, bodyB_file, mesh_file, modes_file,
        sample_file_path, output_file_path,
        max_iters=50
    )
    
except Exception as e:
    print(f"\n Error during execution: {e}")
    import traceback
    traceback.print_exc()
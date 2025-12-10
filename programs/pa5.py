# main for pa5 assignment pa5.py
import os
import sys
from deform_registration import *

# check command line arguments if they were entered ok/are valid
if len(sys.argv) != 2:
    print("Usage: python pa5.py <File_Letter>")
    print("Please entere in the format: python pa5.py A")
    print("Please enter with this format: python pa5.py G")
    sys.exit(1)

# make sure the file lettered entered is a valid argument
file_letter = sys.argv[1].upper()
if file_letter not in "ABCDEFGHJK":
    print(f"You have entered an Invalid file letter '{file_letter}'. Instead please use A-F for debug files or G, H, J, K for the unknown files!.")
    sys.exit(1)
# set up file paths here
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir) 

data_folder = os.path.join(parent_dir, "2025_PA345_Student_Data")
if not os.path.exists(data_folder):
    data_folder = os.path.join(parent_dir, "2025 345 Student Data") 
    
if not os.path.exists(data_folder):
    # exit with an error if the daata folder can't be found
    print(f"the data folder couldn't be found..")
    sys.exit(1)

# make sure file is either debug or unknown
is_debug = file_letter in "ABCDEF"
file_type = "Debug" if is_debug else "Unknown"

sample_filename = f"PA5-{file_letter}-{file_type}-SampleReadingsTest.txt"
output_filename = f"PA5-{file_letter}-{file_type}-Output.txt"

# full file path listed out here
sample_file_path = os.path.join(data_folder, sample_filename)
output_file_path = os.path.join(current_dir, output_filename)

# other input files
bodyA_file = os.path.join(data_folder, "Problem5-BodyA.txt")
bodyB_file = os.path.join(data_folder, "Problem5-BodyB.txt")
mesh_file = os.path.join(data_folder, "Problem5MeshFile.sur")
modes_file = os.path.join(data_folder, "Problem5Modes.txt")

if not os.path.exists(sample_file_path):
    print(f"The Sample file couldn't be found: {sample_file_path}")
    sys.exit(1)

# now call the main function to solve pa5

try:
    s_k, c_k, lambdas = solve_pa5(
        bodyA_file, bodyB_file, mesh_file, modes_file,
        sample_file_path, output_file_path,
        max_iters=50
    )
   #print out results here 
except Exception as e:
    print(f"\n there was an error during the execution of this program: {e}")
    import traceback
    traceback.print_exc()



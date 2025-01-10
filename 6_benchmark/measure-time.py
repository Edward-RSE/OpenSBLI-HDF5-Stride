import sys
from pathlib import Path

if len(sys.argv) != 2:
    print("Please provide a file name to extract times from")
    sys.exit(1)

filename = sys.argv[1]
sim_output = Path(filename).open()
times = {
    "regular": [],
    "slice": [],
    "slice_with_stride": [],
}


for line in sim_output:
    try:
        if "strided" in line:
            times["slice_with_stride"].append(float(line.split()[-1]))
        elif "slice" in line:
            times["slice"].append(float(line.split()[-1]))
        elif ": opensbli_output.h5:" in line:
            times["regular"].append(float(line.split()[-1]))
    except:
        print(line)

sim_output.close()

calculate_average = lambda l: round(sum(l) / len(l) * 1000)
print("Regular          ", calculate_average(times["regular"]), "ms")
print("Slice            ", calculate_average(times["slice"]), "ms")
print("Slice with stride", calculate_average(times["slice_with_stride"]), "ms")


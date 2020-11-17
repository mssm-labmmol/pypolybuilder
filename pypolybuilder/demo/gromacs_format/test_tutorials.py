import subprocess
import os

head = os.getcwd()


sdirs = ["dendrimer/PAMAM",
         "dendrimer/PAMAM_fix_connectors",
         "dendrimer/PPI",
         "dendrimer/PAMAMPPI_Janus_like",
         "dendrimer/SPL7013",
         "polymer/Octane",
         "polymer/PNIPAM",
         "polymer/PNIPAM_fix_connection",
         "polymer/PAMAM_PPI_Half",
         "polymer/PolyEtyleneglycol",
         "polymer/Poly_p-benzamide",
         "dendrimer/test_test_func",
         "non-existent-directory"]


def look_for_gromacs():
    found_gmx = False
    gmx_path = ""
    for bin_name in ["gmx", "gmx_mpi"]:
        session = subprocess.Popen(["which", "{}".format(bin_name)], shell=False, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = session.communicate()
        if stderr:
            continue
        decoded = str(stdout.decode("utf-8"))
        if len(decoded) > 0:
            found_gmx = True
            gmx_path = str(stdout.decode("utf-8")).rstrip("\n")
            break
    return found_gmx, gmx_path


found_gmx, gmx_path = look_for_gromacs()

if found_gmx:
    print("Found gromacs executable:\n{}\n\n".format(gmx_path))
else:
    print("Couldn't find gromacs executable. Will not test gromacs test_functions...\n")

# quit()

print("Testing {} demo directories.\nThis may take a while ...\n".format(len(sdirs)))

outcome = {}
for sdir in sdirs:
    print("##############################")
    print(sdir, "\n", "\n")
    temp_dir = os.path.join(head, sdir)
    if os.path.exists(temp_dir):
        os.chdir(temp_dir)
        # subprocess.call(['sh', 'how_to_run_this_example.txt'])

        session = subprocess.Popen(['sh', 'how_to_run_this_example.txt'], shell=False,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = session.communicate()
        print(str(stdout.decode("utf-8")))

        if stderr:
            print("##############################")
            print("### STDERR ###")
            print("##############################")
            print(sdir, "\n", "\n")
            stderrstring = str(stderr.decode("utf-8"))
            print(stderrstring)
            if "Error" in stderrstring:
                outcome[sdir] = "Exited with error."
                continue
            else:
                outcome[sdir] = "Finished with warning."

        outcome.setdefault(sdir, "Finished successfully")  # if no error entries, it must have been a success...
        # HERE PERFORM GROMACS MINIMIZATION FOR FURTHER TESTING OF TOPOLOGY
    else:
        outcome[sdir] = "Does not exist...?"

print("##############################")
print("##############################")
print("Tested {} sdirs.".format(len(sdirs)))

for sdir, outc in outcome.items():
    print("{:<40}- \t{}".format(sdir, outc))


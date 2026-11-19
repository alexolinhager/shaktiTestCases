#!/usr/bin/env python3
"""Run SHAKTI test cases, then their diagnostics (Python version of runTestCases.m).

Nothing is written in the repository: each test case runs in its own subfolder
of the run directory, RUNDIR/<test case>. The scripts and input files of the
test case are copied there, then runme.py and diagnostics.py are run there.
Results (.nc files) are saved in that subfolder, and videos and figures in its
'diagnostics' subfolder.

From a terminal:
   python runTestCases.py CASES NPROC RUNDIR [RESOLUTION]

   CASES       all, the name of one test case, or names separated by commas:
                  singleMoulin           (about 720 hourly time steps)
                  helheimSeasonal        (3 scenarios of 1 year)
                  thwaitesInterventions  (spin-up, capacity and 3 experiments)
   NPROC       number of processors, used for all the runs
   RUNDIR      run directory, created if needed (outside of the repository)
   RESOLUTION  low (default) or high: mesh of helheimSeasonal and
               thwaitesInterventions (singleMoulin has one mesh only)

Examples:
   python runTestCases.py all 4 ~/shaktiRuns
   python runTestCases.py helheimSeasonal 8 /scratch/shaktiRuns high
   python runTestCases.py singleMoulin,thwaitesInterventions 2 ~/shaktiRuns low

From Python:
   from runTestCases import runTestCases
   runTestCases('all', 4, '~/shaktiRuns')
   runTestCases(['singleMoulin', 'thwaitesInterventions'], 2, '~/shaktiRuns', 'low')
"""
import os
import runpy
import shutil
import sys
import time
import traceback

AVAILABLE = ['singleMoulin', 'helheimSeasonal', 'thwaitesInterventions']

# Input files of each test case, copied to the run directory with runme.py and
# diagnostics.py
INPUTS = {'singleMoulin': ['moulin.py', 'outline.exp'],
          'helheimSeasonal': ['merra2_runoff.mat', 'HelheimWinterSpinup_300m.nc', 'HelheimWinterSpinup_70m.nc'],
          'thwaitesInterventions': ['TIS_island.exp', 'Thwaites_SHAKTI_init.nc', 'Thwaites_SHAKTI_init_500mto10km.nc']}


def runTestCases(cases, nproc, rundir, resolution='low'):
    # Check inputs
    if isinstance(cases, str):
        cases = AVAILABLE if cases.lower() == 'all' else cases.split(',')
    for case in cases:
        if case not in AVAILABLE:
            raise ValueError("runTestCases: test case '" + case + "' not found, choose among: " + ', '.join(AVAILABLE) + ", or 'all'")
    if int(nproc) != nproc or nproc < 1:
        raise ValueError('runTestCases: nproc (number of processors) should be a positive integer')
    if resolution not in ['low', 'high']:
        raise ValueError("runTestCases: resolution should be 'low' or 'high'")
    if not isinstance(rundir, str) or len(rundir) == 0:
        raise ValueError('runTestCases: rundir (run directory) should be the name of a directory')

    # Run directory: absolute path, outside of the repository
    root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    rundir = os.path.realpath(os.path.expanduser(rundir))
    if rundir == root or rundir.startswith(root + os.sep):
        raise ValueError('runTestCases: the run directory should be outside of the repository (' + root + ')')
    os.makedirs(rundir, exist_ok=True)
    if root not in sys.path:
        sys.path.insert(0, root)  # shaktiDiagnostics.py

    # ISSM: load its path if needed
    try:
        import model  # noqa: F401
    except ImportError:
        if os.getenv('ISSM_DIR') is None:
            raise RuntimeError('runTestCases: ISSM is not on the Python path and the ISSM_DIR environment variable is not set')
        sys.path.append(os.getenv('ISSM_DIR') + '/src/m/dev')
        import devpath  # noqa: F401

    # Run the test cases one after the other; if one fails, go on with the next
    origin = os.getcwd()
    status = []
    elapsed = []
    for case in cases:
        folder = os.path.join(rundir, case)
        print('\n=== %s (nproc=%d, resolution=%s) in %s ===' % (case, nproc, resolution, folder))
        start = time.time()
        try:
            os.makedirs(folder, exist_ok=True)
            for file in ['runme.py', 'diagnostics.py'] + INPUTS[case]:
                shutil.copy2(os.path.join(root, case, file), folder)
            # runme.py and diagnostics.py are scripts: they read nproc and resolution from their global variables
            os.chdir(folder)
            runpy.run_path('runme.py', init_globals={'nproc': int(nproc), 'resolution': resolution})
            runpy.run_path('diagnostics.py', init_globals={'resolution': resolution})
            status.append('done')
        except Exception:
            status.append('FAILED')
            traceback.print_exc()
        finally:
            os.chdir(origin)
        elapsed.append(time.time() - start)

    # Summary
    print('\n%-24s %-8s %s' % ('test case', 'status', 'time (min)'))
    for case, s, e in zip(cases, status, elapsed):
        print('%-24s %-8s %.1f' % (case, s, e / 60))
    if 'FAILED' in status:
        raise RuntimeError('runTestCases: %d of %d test cases failed (see above)' % (status.count('FAILED'), len(cases)))


if __name__ == '__main__':
    if len(sys.argv) not in [4, 5]:
        print(__doc__)
        sys.exit('runTestCases: the test cases, the number of processors and the run directory are required')
    runTestCases(sys.argv[1], int(sys.argv[2]), *sys.argv[3:])

#!/usr/bin/env python3
"""Run SHAKTI test cases, then their diagnostics (Python version of runTestCases.m).

For each test case, this runs runme.py and then diagnostics.py in the folder
of the test case. Results (.nc files) are saved in that folder, and videos and
figures in its 'diagnostics' subfolder.

From a terminal:
   python runTestCases.py CASES NPROC [RESOLUTION]

   CASES       all, the name of one test case, or names separated by commas:
                  singleMoulin           (about 720 hourly time steps)
                  helheimSeasonal        (3 scenarios of 1 year)
                  thwaitesInterventions  (spin-up, capacity and 3 experiments)
   NPROC       number of processors, used for all the runs
   RESOLUTION  low (default) or high: mesh of helheimSeasonal and
               thwaitesInterventions (singleMoulin has one mesh only)

Examples:
   python runTestCases.py all 4
   python runTestCases.py helheimSeasonal 8 high
   python runTestCases.py singleMoulin,thwaitesInterventions 2 low

From Python:
   from runTestCases import runTestCases
   runTestCases('all', 4)
   runTestCases(['singleMoulin', 'thwaitesInterventions'], 2, 'low')
"""
import os
import runpy
import sys
import time
import traceback

AVAILABLE = ['singleMoulin', 'helheimSeasonal', 'thwaitesInterventions']


def runTestCases(cases, nproc, resolution='low'):
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

    # ISSM: load its path if needed
    try:
        import model  # noqa: F401
    except ImportError:
        if os.getenv('ISSM_DIR') is None:
            raise RuntimeError('runTestCases: ISSM is not on the Python path and the ISSM_DIR environment variable is not set')
        sys.path.append(os.getenv('ISSM_DIR') + '/src/m/dev')
        import devpath  # noqa: F401

    # Run the test cases one after the other; if one fails, go on with the next
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    origin = os.getcwd()
    status = []
    elapsed = []
    for case in cases:
        print('\n=== %s (nproc=%d, resolution=%s) ===' % (case, nproc, resolution))
        start = time.time()
        try:
            # runme.py and diagnostics.py are scripts: they read nproc and resolution from their global variables
            os.chdir(os.path.join(root, case))
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
    if len(sys.argv) not in [3, 4]:
        print(__doc__)
        sys.exit('runTestCases: the test cases and the number of processors are required')
    runTestCases(sys.argv[1], int(sys.argv[2]), *sys.argv[3:])

from onedep import __apiUrl__
from onedep.api.Validate import Validate
import time
import getopt
import sys
import os

def displayStatus(sD, exitOnError=True):
    if 'onedep_error_flag' in sD and sD['onedep_error_flag']:
        print("OneDep error: %s\n" % sD['onedep_status_text'])
        if exitOnError:
            raise SystemExit()
    else:
        if 'status' in sD:
            print("OneDep status: %s\n" % sD['status'])

def validate_cif(modelFilePath, sfFilePath, prefix):
    val = Validate(apiUrl=__apiUrl__)
    print("initiating new session")
    rD = val.newSession()
    displayStatus(rD)
    print(f"reading input model {modelFilePath}")
    rD = val.inputModelXyzFile(modelFilePath)
    print(f"reading structure factors {sfFilePath}")
    rD = val.inputStructureFactorFile(sfFilePath)
    displayStatus(rD)
    print('running validation; starting...')
    rD = val.run()
    print('status...')
    displayStatus(rD)

    #
    #   Poll for service completion -
    #
    it = 0
    sl = 2
    while (True):
       it += 1
       pause = it * it * sl
       time.sleep(pause)
       rD = val.getStatus()
       if rD['status'] in ['completed', 'failed']:
          print(f" => {rD['status']}")
          break
       print("[%4d] Pausing for %4d (seconds)\n" % (it, pause))

    print("saving pdf and xml reports")
    pdf = f"xray-report-{prefix}.pdf"
    xml = f"xray-report-{prefix}.xml"
    rD = val.getReport(pdf)
    rD = val.getReportData(xml)

    print("finished running validation")


def main(argv):
    prefix = None
    modelFilePath = None
    sfFilePath = None

    try:
        opts, args = getopt.getopt(argv, "m:s:o:h", ["model=", "sf=", "output=", "help"])
    except getopt.GetoptError:
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print("pdb_validate_cif.py -m model.cif -s sf.cif -o output_prefix")
            sys.exit(2)
        elif opt in ("-m", "--model"):
            modelFilePath = os.path.abspath(arg)
        elif opt in ("-s", "--sf"):
            sfFilePath = os.path.abspath(arg)
        elif opt in ("-o", "--output"):
            prefix = arg

    if os.path.isfile(modelFilePath) and os.path.isfile(sfFilePath) and prefix:
        if modelFilePath.endswith('cif') and sfFilePath.endswith('cif'):
            validate_cif(modelFilePath, sfFilePath, prefix)
        else:
            print("error: model and sf files must be in mmcif format")
    else:
        print("pdb_validate_cif.py -m model.cif -s sf.cif -o output_prefix")


if __name__ == "__main__":
    main(sys.argv[1:])

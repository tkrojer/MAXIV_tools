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
    print('status')
    displayStatus(rD)
    print('running validation; starting...')
    rD = val.run()
    print('status')
    displayStatus(rD)

    #
    #   Poll for service completion -
    #
    it = 0
    sl = 2
    while (True):
       #    Pause -
#       print('here')
       it += 1
       pause = it * it * sl
       time.sleep(pause)
       rD = val.getStatus()
       if rD['status'] in ['completed', 'failed']:
          break
       print("[%4d] Pausing for %4d (seconds)\n" % (it, pause))
   #
   #
#    lt = time.strftime("%Y%m%d%H%M%S", time.localtime())
#    fnR = "xray-report-%s.pdf" % lt
##    fnR = "xray-report-%s.xml" % lt
#    rD = val.getReport(fnR)
##    rD = val.getReportData(fnR)

    pdf = f"xray-report-{prefix}.pdf"
    xml = f"xray-report-{prefix}.xml"
    #    fnR = "xray-report-%s.xml" % lt
    rD = val.getReport(pdf)
    rD = val.getReportData(xml)

    print("finished running validation")


def main(argv):
    #    sample_id = None
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
        validate_cif(modelFilePath, sfFilePath, prefix)
    else:
        print("pdb_validate_cif.py -m model.cif -s sf.cif -o output_prefix")


if __name__ == "__main__":
    main(sys.argv[1:])
#    sample_id = 'AR-F2XUniversal-P05D09a_1'
#    modelFilePath = '/Users/tobkro/Scripts/PDB_query/pdb/AR-F2XUniversal-P02H11a_1-pandda-model_edited_refine_001.cif'
#    sfFilePath = '/Users/tobkro/Scripts/PDB_query/pdb/AR-F2XUniversal-P02H11a_1_data_F.cif'


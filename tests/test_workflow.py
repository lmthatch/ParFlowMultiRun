import os
import os.path
import hashlib
from runScripts.SC_GenerateParamSet import main as generate
from runScripts.RunParSet import main as run


def test_output_pfb_files():
    """
    An end-end test that generates 100 distinct parameter sets, runs the (1-indexed) parameter set at index 2,
    and compares the hash of the last pressure/saturation pfb files generated against a previously computed hash.
    This test exists while code refactoring is in progress and will eventually be replaced by more fine-grained tests.
    """
    output_folder = 'test_folder'
    generate([100, 'input_files/SCInputVariables_PFOnly_20200825.csv', output_folder, 10, 123])

    # TODO: The run function changes the current folder on us!
    # Till this is fixed, we explicitly save the current working folder and cd back to it later
    CWD = os.getcwd()

    run([output_folder, 2])

    # TODO: The run function changes the current folder on us!
    # Till this is fixed, we explicitly cd to the folder we were in
    os.chdir(CWD)

    assert hashlib.md5(open(os.path.join(output_folder, 'test1', 'test.out.press.00048.pfb'), 'rb').read())\
               .hexdigest() == 'b059c7442a528950c8d9ecfdfcc5e3d2'

    assert hashlib.md5(open(os.path.join(output_folder, 'test1', 'test.out.satur.00048.pfb'), 'rb').read())\
               .hexdigest() == '05a7f16b1efa53062e155442a4c961b7'

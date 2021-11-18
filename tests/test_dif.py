""" Tests for dif.py """

import os
import platform
import random
import re
import string
from subprocess import getstatusoutput

PRG = './dif.py'
RUN = f'python3 {PRG}' if platform.system() == 'Windows' else PRG
INPUT1 = 'tests/inputs/file1.fa'


# --------------------------------------------------
def test_exists() -> None:
    """ Program exists """

    assert os.path.exists(PRG)


# --------------------------------------------------
def test_usage() -> None:
    """ Prints usage """

    for arg in ['-h', '--help']:
        rv, out = getstatusoutput(f'{RUN} {arg}')
        assert rv == 0
        assert out.lower().startswith('usage:')


# --------------------------------------------------
def test_dies_no_args() -> None:
    """ Dies with no arguments """

    rv, out = getstatusoutput(RUN)
    assert rv != 0
    assert out.lower().startswith('usage:')


# --------------------------------------------------
def test_bad_file() -> None:
    """ Die on bad input file """

    bad = random_filename()
    retval, out = getstatusoutput(f'{RUN} -o foo.png {bad}')
    assert retval != 0
    assert re.match('usage:', out, re.IGNORECASE)
    assert re.search(f"No such file or directory: '{bad}'", out)


# --------------------------------------------------
def test_missing_outfile() -> None:
    """ Die on missing --outfile """

    retval, out = getstatusoutput(f'{RUN} {INPUT1}')
    assert retval != 0
    assert re.match('usage:', out, re.IGNORECASE)
    assert re.search("required: -o/--outfile", out)


# --------------------------------------------------
def test_bad_size() -> None:
    """ Die on missing input """

    for bad_size in [random.randint(0, 999), random.randint(30001, 50000)]:
        retval, out = getstatusoutput(
            f'{RUN} -o foo.png --size {bad_size} {INPUT1}')
        assert retval != 0
        assert re.match('usage:', out, re.IGNORECASE)
        assert re.search(fr"error: --size \({bad_size}\) must be between", out)


# --------------------------------------------------
def random_filename() -> str:
    """ Generate a random filename """

    return ''.join(random.choices(string.ascii_uppercase + string.digits, k=5))

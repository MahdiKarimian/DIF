#!/usr/bin/env python3
"""
Author : Mahdi Karimian
Purpose: Convert DNA sequence to image
"""

from Bio import SeqIO
# from matplotlib.pyplot import imshow
import numpy as np
import cv2  # python-opencv
# import random
import argparse
from typing import Dict, NamedTuple, TextIO, Tuple

# pylint: disable=no-member


class Args(NamedTuple):
    """ Command-line arguments """
    file: TextIO
    outfile: str
    size: int
    recordid: str
    plotsize: int


class Options(NamedTuple):
    """ Settings """
    image_size: int  # Output image size, 30000 for big sequence
    gain_xy: float  # Zoom the sequence, this parameter is in development phase
    seq_color: dict[str, Tuple[int, int, int, int]]  # Color for A, C, T, G
    max_nucl_dist: int  # Advance param: max distance of two nucletide in graph


# class Steps(Dict[str, Tuple[int, int]]):
#     """ Definition of how coords should change by base """


# --------------------------------------------------
def get_args() -> Args:
    """Get command-line arguments"""
    parser = argparse.ArgumentParser(
        description='Convert DNA sequence to image')

    parser.add_argument('file',
                        type=argparse.FileType('rt'),
                        metavar='FILE',
                        help='FASTA input file')

    parser.add_argument('-o',
                        '--outfile',
                        type=str,
                        metavar='FILE',
                        help='Image output file name',
                        required=True)

    parser.add_argument('-s',
                        '--size',
                        type=int,
                        metavar='INT',
                        help='Image output size (1000-30000)',
                        default=1000)

    parser.add_argument('-r',
                        '--recordid',
                        type=str,
                        metavar='SEQ',
                        default='*',
                        help='Record ID')

    parser.add_argument('-p',
                        '--plotsize',
                        type=int,
                        metavar='INT',
                        default=10000000,
                        help='Limit the size of plot default 10M')

    args = parser.parse_args()

    if not 1000 <= args.size <= 30000:
        parser.error(f'--size ({args.size}) must be between 1000 and 30,000.')

    return Args(args.file, args.outfile, args.size, args.recordid,
                args.plotsize)


# --------------------------------------------------
def main() -> None:
    """The Good Stuff"""

    args = get_args()

    option = Options(
        args.size, 0.2, {
            'a': (255, 0, 0, 255),
            'c': (0, 255, 0, 255),
            't': (0, 0, 255, 255),
            'g': (128, 128, 128, 255)
        }, 5)

    # if you want another seq
    for record in SeqIO.parse(args.file, "fasta"):
        if record.id == args.recordid:
            print(record.id)
            break
        if args.recordid == '*':
            break

    # record.seq
    len(record.seq)
    print(f"Select this record {record.id} ")

    # generate image footprint for seq and save as test.png file

    image = generate_footprint(record.seq[:args.plotsize], option)
    cv2.imwrite(str(args.outfile), image)
    print("end of Generating image")


# --------------------------------------------------
def generate_footprint(seq: str, option: Options):
    """ Generate image from sequence """

    init_image_size = option.image_size
    seq_color = option.seq_color
    inc_val = option.max_nucl_dist
    gain = option.gain_xy
    inc_val = int(gain * inc_val)

    max_x = init_image_size
    min_x = 0

    max_y = init_image_size
    min_y = 0

    width, height = init_image_size, init_image_size
    x1, y1 = int(init_image_size / 2), int(init_image_size / 2)
    image = np.zeros((height, width, 4), np.uint8)

    # Variables to find many same nucleotide in series TODO
    # now deleted after 5000 sample
    prev_base = ''
    st_count = 0

    steps = get_steps(inc_val)

    for base in seq.lower():

        x2, y2 = get_coords(steps, base, (x1, y1))

        if base == prev_base:
            st_count += 1
        else:
            st_count = 0

        if st_count == 100:
            print(f'st_count {base} reached 100.')
            continue

        prev_base = base

        # dont_plot = 0

        if max_x < x2:
            max_x = x2
            x2 = 0 + 50
            # dontplot = 1

        if max_y < y2:
            max_y = y2
            y2 = 0 + 50
            # dontplot = 1

        if min_x > x2:
            min_x = x2
            x2 = max_x - 50
            # dontplot = 1

        if min_y > y2:
            min_y = y2
            y2 = max_y - 50
            # dontplot = 1

        line_thickness = 3

        if np.abs(x2 - x1) < 100 and np.abs(y2 - y1) < 100:
            image = cv2.line(image, (x1, y1), (x2, y2),
                             color=seq_color[base],
                             thickness=line_thickness)
        x1, y1 = x2, y2

    print("x min max y min max", min_x, ',', max_x, ',', min_y, ',', max_y)
    return image
    # imshow(image)


# --------------------------------------------------
def get_steps(inc: int) -> Dict[str, Tuple[int, int]]:
    """ Generate step dictionary from increment"""

    # Steps dictionary where 'nuc': (dx, dy)
    steps = {
        'a': (inc * 2, inc),
        'c': (-inc, inc),
        't': (-inc * 2, -inc),
        'g': (inc, -inc),
    }

    return steps


# --------------------------------------------------
def test_get_steps() -> None:
    """ Test get_steps """

    assert get_steps(1) == {
        'a': (2, 1),
        'c': (-1, 1),
        't': (-2, -1),
        'g': (1, -1)
    }
    assert get_steps(2) == {
        'a': (4, 2),
        'c': (-2, 2),
        't': (-4, -2),
        'g': (2, -2)
    }


# --------------------------------------------------
def get_coords(steps: Dict[str, Tuple[int, int]], base: str,
               current: Tuple[int, int]) -> Tuple[int, int]:
    """ Get new coordinates based on nucleotide """

    x1, y1 = current

    dx, dy = steps.get(base, (0, 0))

    return x1 + dx, y1 + dy


# --------------------------------------------------
def test_get_coords() -> None:
    """ Test get_coords """

    steps = get_steps(1)
    
    # Unrecognized residue
    assert get_coords(steps, '', (0, 0)) == (0, 0)
    assert get_coords(steps, 'm', (0, 0)) == (0, 0)
    assert get_coords(steps, '', (5, 5)) == (5, 5)

    # Standard residues, increment = 1
    assert get_coords(steps, 'a', (0, 0)) == (2, 1)
    assert get_coords(steps, 'c', (0, 0)) == (-1, 1)
    assert get_coords(steps, 't', (0, 0)) == (-2, -1)
    assert get_coords(steps, 'g', (0, 0)) == (1, -1)

    # Standard residues, increment = 2
    steps = get_steps(2)
    assert get_coords(steps, 'a', (0, 0)) == (4, 2)
    assert get_coords(steps, 'c', (0, 0)) == (-2, 2)
    assert get_coords(steps, 't', (0, 0)) == (-4, -2)
    assert get_coords(steps, 'g', (0, 0)) == (2, -2)


# --------------------------------------------------
if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Author : Mahdi Karimian
Purpose: Convert DNA sequence to image
"""

from Bio import SeqIO
from matplotlib.pyplot import imshow
import numpy as np
import cv2  # python-opencv
import random
import argparse
from typing import List, NamedTuple, TextIO


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
    seq_color: List[tuple]  # Color for 'A' , 'C', 'T', 'G'
    max_nucl_distance_px: int  # Advance parameter: max distance of two nucletide in graph


# --------------------------------------------------
def get_args():
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
def main():
    """The Good Stuff"""

    args = get_args()

    option = Options(args.size, 0.2, [(255, 0, 0, 255), (0, 255, 0, 255),
                                      (0, 0, 255, 255), (128, 128, 128, 255)],
                     5)

    # if you want another seq
    for record in SeqIO.parse(args.file, "fasta"):
        if (record.id == args.recordid):
            print(record.id)
            break
        elif (args.recordid == '*'):
            break

    #record.seq
    len(record.seq)
    print("Select this record {} ".format(record.id))

    #generate image footprint for seq and save as test.png file

    image = generate_footprint(record.seq[:args.plotsize], option)
    cv2.imwrite(str(args.outfile), image)
    print("end of Generating image")


# --------------------------------------------------
def generate_footprint(seq, option):
    """ Generate image from sequence """

    init_image_size = option.image_size
    seq_color = option.seq_color
    inc_val = option.max_nucl_distance_px
    gain = option.gain_xy
    inc_val = int(gain * inc_val)

    max_x = init_image_size
    min_x = 0

    max_y = init_image_size
    min_y = 0

    width, height = init_image_size, init_image_size
    x1, y1 = int(init_image_size / 2), int(init_image_size / 2)
    image = np.zeros((height, width, 4), np.uint8)

    state = 'a'  #satte for many same nucletide in series TODO now deleted after 5000 sample
    st_count = 0

    for c in seq:

        if (c == 'a' or c == 'A'):
            in_val = 0
            if (state == 'a'):
                st_count = st_count + 1
            else:
                st_count = 0
            if (st_count == 100):
                print("st_count a reach 5000")
                continue
        elif (c == 'c' or c == 'C'):
            in_val = 1
            if (state == 'c'):
                st_count = st_count + 1
            else:
                st_count = 0
            if (st_count == 100):
                print("st_count c reach 5000")
                continue
        elif (c == 't' or c == 'T'):
            in_val = 2
            if (state == 't'):
                st_count = st_count + 1
            else:
                st_count = 0
            if (st_count == 100):
                print("st_count t reach 5000")
                continue

        elif (c == 'g' or c == 'G'):
            in_val = 3
            if (state == 'g'):
                st_count = st_count + 1
            else:
                st_count = 0
            if (st_count == 100):
                print("st_count g reach 5000")
                continue
        else:
            #print('error char is ',c)
            continue

        #print(in_val)
        if (in_val == 0):
            x2, y2 = x1 + int(inc_val * 2), y1 + inc_val
        elif (in_val == 1):
            x2, y2 = x1 - inc_val, y1 + inc_val
        elif (in_val == 2):
            x2, y2 = x1 - int(inc_val * 2), y1 - inc_val
        elif (in_val == 3):
            x2, y2 = x1 + inc_val, y1 - inc_val

        dont_plot = 0

        if (max_x < x2):
            max_x = x2
            x2 = 0 + 50
            #dontplot = 1

        if (max_y < y2):
            max_y = y2
            y2 = 0 + 50
            #dontplot = 1

        if (min_x > x2):
            min_x = x2
            x2 = max_x - 50
            #dontplot = 1

        if (min_y > y2):
            min_y = y2
            y2 = max_y - 50
            #dontplot = 1

        line_thickness = 3

        if (np.abs(x2 - x1) < 100 and np.abs(y2 - y1) < 100):
            image = cv2.line(image, (x1, y1), (x2, y2),
                             color=seq_color[in_val],
                             thickness=line_thickness)
        x1, y1 = x2, y2

    print("x min max y min max", min_x, ',', max_x, ',', min_y, ',', max_y)
    return image
    #imshow(image)


# --------------------------------------------------
if __name__ == '__main__':
    main()
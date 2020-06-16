# All code of gfapy is released under the following ISC license.
# It is functionally equivalent to a two-term BSD copyright with
# language removed that is made unnecessary by the Berne convention.
# See http://openbsd.org/policy.html for more information on copyrights.

# Copyright (c) 2017 Giorgio Gonnella and CONTRIBUTORS
# Copyright (c) 2017 Center for Bioinformatics, University of Hamburg

# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.

# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

import argparse
import sys
import random

op = argparse.ArgumentParser(description=__doc__)
op.add_argument("--segments", "-s", type=int,
                help="number of segments", required=True)
op.add_argument("--slen", "-l", type=int, default=100,
                help="lenght of segments sequence")
op.add_argument("--with-sequence", "-w", action="store_true")
op.add_argument("--dovetails-per-segment", "-d",
                help="average number of dovetail edges per segment",
                default=2.0, type=float)
op.add_argument('--gfa-version', "-g", default="gfa1",
                help="gfa version", choices=("gfa1", "gfa2"))
op.add_argument('--version', action='version', version='%(prog)s 1.0')
opts = op.parse_args()

if opts.segments < 0:
    sys.stderr.write("Error: the number of segments must be "+
                    ">= 0 ({})\n".format(opts.segments))
    exit(1)
if opts.dovetails_per_segment < 0:
    sys.stderr.write("Error: the average number of dovetails per segment must "+
                    "be >= 0 ({})\n".format(opts.dovetails_per_segment))
    exit(1)
if opts.slen <= 0:
    sys.stderr.write("Error: the length of segments sequence must be > 0"+
                   " ({})\n".format(opts.slen))
    exit(1)

if opts.gfa_version == "gfa1":
    print("H\tVN:Z:1.0")
else:
    print("H\tVN:Z:2.0")

def random_sequence(slen):
    sequence = []
    for i in range(slen):
        sequence.append(random.choice('ACGT'))
    return "".join(sequence)

for i in range(opts.segments):
    if opts.with_sequence:
        sequence = random_sequence(opts.slen)
    else:
        sequence = "*"
    if opts.gfa_version == "gfa1":
        print("S\ts{}\t{}\tLN:i:{}".format(i, sequence, opts.slen))
    else:
        print("S\ts{}\t{}\t{}".format(i, opts.slen, sequence))

n_dovetails = int(opts.segments * opts.dovetails_per_segment)
edges = {}
for i in range(n_dovetails):
    edge = False
    while not edge:
        s_from = random.randint(0, opts.segments-1)
        s_from_or = random.choice('+-')
        s_to = random.randint(0, opts.segments-1)
        s_to_or = random.choice('+-')
        if s_from not in edges:
            edges[s_from] = {'+': {}, '-': {}}
        if s_to not in edges[s_from][s_from_or]:
            edges[s_from][s_from_or][s_to] = {'+': False, '-': False}
        if not edges[s_from][s_from_or][s_to][s_to_or]:
            edges[s_from][s_from_or][s_to][s_to_or] = True
            edge = True
    ovlen = opts.slen//10
    if ovlen == 0: ovlen = 1
    cigar = "{}M".format(ovlen)
    if opts.gfa_version == "gfa1":
        print("L\ts{}\t{}\ts{}\t{}\t{}\tID:Z:e{}".format(s_from, s_from_or, s_to,
                                                   s_to_or, cigar, i))
    else:
        s_from_begin = opts.slen - ovlen if s_from_or == "+" else 0
        s_from_end = "{}$".format(opts.slen) if s_from_or == "+" else ovlen
        s_to_begin = opts.slen - ovlen if s_to_or == "-" else 0
        s_to_end = "{}$".format(opts.slen) if s_to_or == "-" else ovlen
        print("E\te{}\ts{}{}\ts{}{}\t{}\t{}\t{}\t{}\t{}".format(
            i, s_from, s_from_or, s_to, s_to_or, s_from_begin, s_from_end,
            s_to_begin, s_to_end, cigar))

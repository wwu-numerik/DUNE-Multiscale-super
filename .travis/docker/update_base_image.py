#!/usr/bin/env python3

import os
import subprocess
import sys


def update(branch):
    thisdir = os.path.dirname(os.path.abspath(__file__))
    gdt_super_dir = os.path.join(thisdir, '..', '..',)
    os.chdir(gdt_super_dir)

    subprocess.check_call(['docker', 'build', '-f', '.travis/docker/testing-base/Dockerfile',
                        '-t', 'dunecommunity/testing-base:multiscale_pre_2.4', '.'])
    subprocess.check_call(['docker', 'push', 'dunecommunity/testing-base:multiscale_pre_2.4'])

if __name__ == '__main__':
    if len(sys.argv) > 1:
        branch = sys.argv[1]
    else:
        branch = 'master'

    update(branch)



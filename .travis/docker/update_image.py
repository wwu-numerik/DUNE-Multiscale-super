#!/usr/bin/env python3

import os
import subprocess
import sys

cc_mapping = {'gcc': 'g++', 'clang': 'clang++'}

def update(branch, cc):
    gdt_super_dir = os.path.join(os.path.dirname(__file__), '..', '..',)
    os.chdir(gdt_super_dir)
    cxx = cc_mapping[cc]

    #for sub in [ 'dune-multiscale-dev']:
    for sub in ['dune-multiscale-testing', 'dune-multiscale-dev']:
        dockerfile = os.path.join(os.path.dirname(__file__), sub, 'Dockerfile')
        repo = 'dunecommunity/{}_{}'.format(sub, cc)
        for tag in [branch, 'latest']:
            subprocess.check_call(['docker', 'build', '-f', dockerfile,
                                '-t', '{}:{}'.format(repo, tag), '--build-arg', 'cc={}'.format(cc),
                                '--build-arg', 'cxx={}'.format(cxx), '.'])
        subprocess.check_call(['docker', '--log-level="debug"', 'push', repo])

if __name__ == '__main__':
    if len(sys.argv) > 2:
        ccs = [sys.argv[1]]
        branches = [sys.argv[2]]
    else:
        ccs = list(cc_mapping.keys())
        branches = ['master']

    for b in branches:
        for c in ccs:
            update(b, c)

#!/usr/bin/env python3

"""
Update docker images and templated scripts in dune-xt-*

Usage:
  update.py [options]

Options:
   -h --help       Show this message.
   -v --verbose    Set logging level to debug.
"""

from docopt import docopt
from os import path
import importlib
import os
import stat
import contextlib
from string import Template as stringTemplate
import subprocess
import sys
import logging
import tempfile
import time
import six
import re
try:
    import docker
except ImportError:
    print('missing module: pip install docker')
    sys.exit(1)
from docker.utils.json_stream import json_stream

TAG_MATRIX = {'debian_gcc_full': {'cc': 'gcc', 'cxx': 'g++', 'deletes': "", 'base': 'debian'},
        'debian_clang_full': {'cc': 'clang', 'cxx': 'clang++', 'deletes':"", 'base': 'debian'},
        'arch_gcc_full': {'cc': 'gcc', 'cxx': 'g++', 'deletes': "", 'base': 'arch'},}


@contextlib.contextmanager
def remember_cwd(dirname):
    curdir = os.getcwd()
    try:
        os.chdir(dirname)
        yield curdir
    finally:
        os.chdir(curdir)


class Timer(object):
    def __init__(self, section, log):
        self._section = section
        self._start = 0
        self._log = log
        self.time_func = time.time

    def start(self):
        self.dt = -1
        self._start = self.time_func()

    def stop(self):
        self.dt = self.time_func() - self._start

    def __enter__(self):
        self.start()

    def __exit__(self, type_, value, traceback):
        self.stop()
        self._log('Execution of {} took {} (s)'.format(self._section, self.dt))


class CommitMessageMissing(RuntimeError): pass


def _docker_build(client, **kwargs):
    resp = client.api.build(**kwargs)
    if isinstance(resp, six.string_types):
        return client.images.get(resp)
    last_event = None
    image_id = None
    output = []
    for chunk in json_stream(resp):
        if 'error' in chunk:
            msg = chunk['error'] + '\n' + ''.join(output)
            raise docker.errors.BuildError(msg)
        if 'stream' in chunk:
            output.append(chunk['stream'])
            match = re.search(
                r'(^Successfully built |sha256:)([0-9a-f]+)$',
                chunk['stream']
            )
            if match:
                image_id = match.group(2)
        last_event = chunk
    if image_id:
        return client.images.get(image_id)
    raise docker.errors.BuildError(last_event or 'Unknown')


def _cmd(cmd, logger):
    logger.debug(' '.join(cmd))
    try:
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, universal_newlines=True)
        logger.debug(out)
    except subprocess.CalledProcessError as cp:
        logger.error(cp.output)
        logger.error('Failed: {}'.format(' '.join(cmd)))
        logger.error('Make sure the pushers group has write access to this repo on hub.cocker.com!')
        raise cp


def _build_base(scriptdir, distro, cc, cxx, commit, refname):
    client = docker.from_env(version='auto')
    slug_postfix = '{}_{}'.format(distro, cc)
    logger = logging.getLogger('{}'.format(slug_postfix))
    dockerdir = path.join(scriptdir, 'dune-multiscale-testing')
    dockerfile = path.join(dockerdir, 'Dockerfile')
    repo = 'dunecommunity/dune-multiscale-testing_{}'.format(slug_postfix)
    tag = '{}:{}'.format(repo, commit)
    logger.error(tag)
    with Timer('docker build ', logger.info):
        buildargs = {'COMMIT': commit, 'CC': cc, 'CXX': cxx, 'BASE': distro}
        img = _docker_build(client, rm=False, buildargs=buildargs,
                            tag=tag, path=dockerdir)
        img.tag(repo, refname)
    with Timer('docker push ', logger.info):
        client.images.push(repo)
    return img



if __name__ == '__main__':
    arguments = docopt(__doc__)
    level = logging.DEBUG if arguments['--verbose'] else logging.INFO
    logging.basicConfig(level=level)
    scriptdir = path.dirname(path.abspath(__file__))
    superdir = path.join(scriptdir, '..', '..')

    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], universal_newlines=True).strip()
    commit = os.environ.get('CI_COMMIT_SHA', head)
    refname = os.environ.get('CI_COMMIT_REF_NAME', 'master').replace('/', '_')

    all_compilers = {(f['base'], f['cc'], f['cxx']) for f in TAG_MATRIX.values()}
    base_imgs = [_build_base(scriptdir, base, cc, cxx, commit, refname) for base, cc, cxx in all_compilers]
    client = docker.from_env(version='auto')

